# Model setup: the grid from the input data, its spacing and orientation (AeoLiS
# `gridparams`), and the run-time values AeoLiS keeps in its parameter dictionary.

"""
    AeolisModel{T}

Everything a time step needs besides the state: the parameters, the grid dimensions, the
forcing series, masks, caches for shear computations, the random number generator and
counters. Built by [`AeolisProblem`](@ref)'s integrator; users normally reach it as
`integ.model`.
"""
mutable struct AeolisModel{T<:AbstractFloat,R<:Random.AbstractRNG}
    p::AeolisParameters
    nx::Int
    ny::Int
    nl::Int
    nf::Int
    alpha::Float64               # grid orientation with respect to the x-axis [deg]
    xgrid::Matrix{T}
    ygrid::Matrix{T}
    bed::Matrix{T}
    ne::Matrix{T}
    grain_size::Vector{T}
    grain_dist::Matrix{T}        # one row per layer, or a single row
    wind::Union{Nothing,Matrix{T}}
    wind_speed::Any
    wind_sin::Any
    wind_cos::Any
    tide::Any
    wave_height::Any
    wave_period::Any
    meteo::Any
    masks::Dict{Symbol,Any}
    supply::Union{Nothing,Matrix{T}}
    bedcomp::Union{Nothing,Array{T,4}}
    shear::Any                   # WindShear or RotationGrid for 2-D shear
    okin::Any
    meteo_values::Dict{Symbol,Float64}
    z0_init::Any
    dt_opt::Float64
    dt_prev::Float64
    time::Float64
    rng::R
    counters::Dict{Symbol,Int}
    prev::AeolisState{T}          # the state at the start of the step (AeoLiS `l`)
end

function _count!(m::AeolisModel, name::Symbol, n::Int=1)
    return (m.counters[name] = get(m.counters, name, 0) + n)
end

_asfloat(::Type{T}, a::AbstractArray) where {T} = T.(a)

# Turn an input array into a grid of size (nx+1, ny+1): matrices are given x first; vectors
# are transects along x.
function _as_grid(::Type{T}, a, nxp::Int, nyp::Int, name::Symbol) where {T}
    a === nothing && return fill(T(NaN), nxp, nyp)
    if a isa Number
        return fill(T(a), nxp, nyp)
    elseif ndims(a) == 1
        length(a) == nxp ||
            throw(ConfigurationError("$name has $(length(a)) values; the grid has $nxp points in x"))
        return repeat(T.(a), 1, nyp)
    else
        size(a) == (nxp, nyp) ||
            throw(ConfigurationError("$name has size $(size(a)); the grid is ($nxp, $nyp)"))
        return T.(a)
    end
end

function _check_inputs(p::AeolisParameters)
    p.xgrid_file === nothing && throw(ConfigurationError("no grid: xgrid_file is not set"))
    p.bed_file === nothing && throw(ConfigurationError("no bed: bed_file is not set"))
    if p.ygrid_file === nothing && ndims(p.xgrid_file) == 1
        # AeoLiS requires a y-grid for transects too; a transect at y = 0 is the natural default.
    elseif p.ygrid_file === nothing
        throw(ConfigurationError("a two-dimensional grid needs ygrid_file"))
    end
    checks = [(:wind_file, 3), (:tide_file, 2), (:meteo_file, 6)]
    (p.process_wave || p.process_dune_erosion) && push!(checks, (:wave_file, 3))
    for (name, ncol) in checks
        a = getfield(p, name)
        a === nothing && continue
        (ndims(a) == 2 && size(a, 2) >= ncol) ||
            throw(ConfigurationError("$name must be a matrix with at least $ncol columns, got size $(size(a))"))
    end
    if p.solver in (:steadystate, :steadystatepieter) && length(p.grain_size) > 1
        throw(ConfigurationError("the steadystate solver is not compatible with multiple sediment fractions; use a single fraction or the trunk or pieter solver"))
    end
    _check_option(:scheme, p.scheme, (:euler_forward, :euler_backward, :crank_nicolson))
    _check_option(:solver, p.solver, (:trunk, :pieter, :steadystate, :steadystatepieter))
    _check_option(:wind_convention, p.wind_convention, (:nautical, :cartesian))
    _check_option(:method_grainspeed, p.method_grainspeed,
                  (:windspeed, :duran, :duran_full, :constant))
    _check_option(:method_transport, p.method_transport,
                  (:bagnold, :bagnold_gs, :kawamura, :lettau, :dk, :sauermann,
                   :vanrijn_strypsteen))
    _check_option(:method_roughness, p.method_roughness,
                  (:constant, :constant_nikuradse, :mean_grainsize_initial,
                   :mean_grainsize_adaptive, :median_grainsize_adaptive,
                   :vanrijn_strypsteen))
    _check_option(:method_shear, p.method_shear, (:fft, :quasi2d, :duna2d))
    _check_option(:vegshear_type, p.vegshear_type, (:raupach, :okin))
    _check_option(:veg_model, p.veg_model, (:legacy, :marram))
    p.process_vegetation && p.veg_model == :marram &&
        throw(UnsupportedOptionError(:veg_model, :marram,
                                     "the marram grass grid model is not implemented yet; the single-cell model is (simulate_marram)"))
    _check_option(:boundary_offshore, p.boundary_offshore,
                  (:flux, :noflux, :constant, :uniform, :gradient, :circular, :re_circular))
    _check_option(:boundary_onshore, p.boundary_onshore,
                  (:flux, :noflux, :constant, :uniform, :gradient, :circular))
    _check_option(:boundary_lateral, p.boundary_lateral,
                  (:circular, :constant, :uniform, :gradient, :flux, :re_circular))
    p.th_drylayer && throw(UnsupportedOptionError(:th_drylayer, true,
                                                  "the dry-layer threshold is not implemented in AeoLiS"))
    p.th_humidity && throw(UnsupportedOptionError(:th_humidity, true,
                                                  "the humidity threshold is not implemented in AeoLiS"))
    return nothing
end

function _grain_dist_matrix(::Type{T}, gd) where {T}
    return gd isa AbstractVector ? reshape(T.(gd), 1, :) : T.(gd)
end

"""
    setup_model(p::AeolisParameters; T=Float64, rng=Random.Xoshiro(0))

The model and its initial state for the parameters `p`: the grid (a transect becomes a
quasi-2-D grid of three identical rows, as in AeoLiS), grid spacing and orientation, bed
composition, wind, vegetation and fences. Returns `(model, state)`.
"""
function setup_model(p::AeolisParameters; T::Type{<:AbstractFloat}=Float64,
                     rng::Random.AbstractRNG=Random.Xoshiro(0))
    _check_inputs(p)
    xg = p.xgrid_file
    if ndims(xg) == 1
        # quasi 2-D: three identical rows, y shifted by the x-spacing (AeoLiS `initialize`)
        qnr = 3
        nxp = length(xg)
        xgrid = repeat(T.(xg), 1, qnr)
        dy = xgrid[3, 2] - xgrid[2, 2]
        y0 = p.ygrid_file === nothing ? zeros(T, nxp) : T.(vec(p.ygrid_file))
        ygrid = hcat((y0 .+ i * dy for i in 0:(qnr - 1))...)
        length(p.bed_file) == nxp ||
            throw(ConfigurationError("bed_file has $(length(p.bed_file)) values; the grid has $nxp points"))
        bed = repeat(T.(vec(p.bed_file)), 1, qnr)
    else
        nxp = size(xg, 1)
        xgrid = T.(xg)
        ygrid = T.(p.ygrid_file)
        bed = _as_grid(T, p.bed_file, nxp, size(xg, 2), :bed_file)
    end
    nyp = size(xgrid, 2)
    nx, ny = nxp - 1, nyp - 1
    nl = p.nlayers
    nf = length(p.grain_size)
    ne = _as_grid(T, p.ne_file, nxp, nyp, :ne_file)
    s = AeolisState{T}(nx, ny, nl, nf)
    masks = Dict{Symbol,Any}()
    for name in
        (:wave_mask, :tide_mask, :runup_mask, :threshold_mask, :gw_mask, :vver_mask)
        v = getfield(p, name)
        masks[name] = v === nothing ? 1.0 :
                      (ndims(v) == 1 ? repeat(vec(v), 1, nyp) : reshape(v, nxp, nyp))
    end
    supply = p.supply_file === nothing ? nothing :
             _as_grid(T, p.supply_file, nxp, nyp, :supply_file)
    bedcomp = p.bedcomp_file === nothing ? nothing :
              _bedcomp_array(T, p.bedcomp_file, nxp, nyp, nl, nf)
    model = AeolisModel{T,typeof(rng)}(p, nx, ny, nl, nf, 0.0, xgrid, ygrid, bed, ne,
                                       T.(p.grain_size),
                                       _grain_dist_matrix(T, p.grain_dist),
                                       nothing, nothing, nothing, nothing, nothing, nothing,
                                       nothing, nothing, masks, supply, bedcomp, nothing,
                                       nothing,
                                       Dict{Symbol,Float64}(), nothing, 0.0, 0.0, p.tstart,
                                       rng,
                                       Dict{Symbol,Int}(), copy(s))
    _setup_forcing!(model)
    initialize_grid!(s, model)
    initialize_bed!(s, model)
    initialize_wind!(s, model)
    initialize_vegetation!(s, model)
    initialize_fences!(s, model)
    copyto_state!(model.prev, s)
    return model, s
end

# AeoLiS reshapes bedcomp_file to (ny+1, nx+1, nlayers, nfractions), or to
# (nx+1, nlayers, nfractions) repeated over three rows for a transect.
function _bedcomp_array(::Type{T}, a, nxp, nyp, nl, nf) where {T}
    out = zeros(T, nxp, nyp, nl, nf)
    ndims(a) == 4 && return copyto!(out, T.(a))
    # read_config transposes spatial files (text rows are y), so `vec` is the reading order
    v = vec(a)
    if nyp == 3 && length(v) == nxp * nl * nf
        r = permutedims(reshape(v, nf, nl, nxp), (3, 2, 1))
        for j in 1:3
            out[:, j, :, :] .= r
        end
    else
        length(v) == nxp * nyp * nl * nf ||
            throw(ConfigurationError("bedcomp_file has $(length(v)) values; expected $(nxp*nyp*nl*nf)"))
        r = reshape(v, nf, nl, nxp, nyp)
        out .= permutedims(r, (3, 4, 2, 1))
    end
    return out
end

function _setup_forcing!(m::AeolisModel{T}) where {T}
    p = m.p
    if p.wind_file !== nothing
        w = T.(p.wind_file)
        if p.wind_convention == :nautical
            w[w[:, 3] .== 0, 3] .= 0.01
        else
            w[w[:, 3] .== 270, 3] .= 270.01
            w[:, 3] .= 270.0 .- w[:, 3]
        end
        m.wind = w
        d = w[:, 3] ./ 180.0 .* pi
        m.wind_speed = CircularStepSeries(w[:, 1], w[:, 2])
        m.wind_sin = CircularStepSeries(w[:, 1], sin.(d))
        m.wind_cos = CircularStepSeries(w[:, 1], cos.(d))
    end
    if p.tide_file !== nothing
        m.tide = CircularSeries(p.tide_file[:, 1], p.tide_file[:, 2])
    end
    if p.wave_file !== nothing && ndims(p.wave_file) == 2 && size(p.wave_file, 2) >= 3
        m.wave_height = CircularSeries(p.wave_file[:, 1], p.wave_file[:, 2])
        m.wave_period = CircularSeries(p.wave_file[:, 1], p.wave_file[:, 3])
    end
    if p.meteo_file !== nothing
        m.meteo = [CircularSeries(p.meteo_file[:, 1], p.meteo_file[:, k])
                   for k in 2:size(p.meteo_file, 2)]
    end
    return m
end

"""
    initialize_grid!(s, m)

Coordinates, cell sizes and areas, and the grid orientation `m.alpha` (AeoLiS
`gridparams.initialize`). As in AeoLiS, `dn` is the distance between neighbouring points
along x and `ds` along y; the two are equal on square grids and on transects.
"""
function initialize_grid!(s::AeolisState, m::AeolisModel)
    s.x .= m.xgrid
    s.y .= m.ygrid
    dx = s.x[2, 1] - s.x[1, 1]
    dy = s.y[2, 1] - s.y[1, 1]
    alpha = dx == 0 ? 90.0 : rad2deg(atan(dy / dx))
    (dx <= 0 && dy <= 0) && (alpha += 180.0)
    m.alpha = alpha
    xr, yr = rotate(s.x, s.y, alpha; origin=(mean(s.x), mean(s.y)))
    s.dn .= sqrt((yr[2, 1] - yr[1, 1])^2 + (xr[2, 1] - xr[1, 1])^2)
    s.ds .= sqrt((xr[1, 2] - xr[1, 1])^2 + (yr[1, 2] - yr[1, 1])^2)
    s.dsdn .= s.ds .* s.dn
    s.dsdni .= 1 ./ s.dsdn
    if (s.x[2, 1] - s.x[1, 1]) <= 0 || (s.y[1, 2] - s.y[1, 1]) <= 0
        @warn "after rotation to the horizontal both x and y should be ascending; the solver may give wrong results"
    end
    return s
end

"""
    grid_rotate!(s, m, angle)

Rotate the coordinates and every vector field (shear stress, shear velocity, wind, fluxes,
grain velocities) by `angle` degrees, and add `angle` to the wind direction (AeoLiS
`grid_rotate`).
"""
function grid_rotate!(s::AeolisState, m::AeolisModel, angle::Real)
    x, y = rotate(s.x, s.y, angle; origin=(mean(s.x), mean(s.y)))
    s.x .= x
    s.y .= y
    for (a, b) in ((:taus, :taun), (:taus0, :taun0), (:ustars, :ustarn),
                   (:ustars0, :ustarn0), (:uws, :uwn))
        u, v = rotate(getfield(s, a), getfield(s, b), angle)
        getfield(s, a) .= u
        getfield(s, b) .= v
    end
    for i in 1:m.nf
        for (a, b) in ((:qs, :qn), (:us, :un))
            A, B = getfield(s, a), getfield(s, b)
            u, v = rotate(A[:, :, i], B[:, :, i], angle)
            A[:, :, i] .= u
            B[:, :, i] .= v
        end
    end
    s.udir .+= angle
    return s
end
