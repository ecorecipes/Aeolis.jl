# Model state. One table lists every spatial variable with its dimensions, unit and
# description; `AeolisState` is generated from it. Arrays are (nx+1, ny+1[, nlayers][,
# nfractions]): x varies fastest, so `vec` of a grid has the cell order of AeoLiS.

const XY = (:nx, :ny)
const XYF = (:nx, :ny, :nfractions)
const XYL = (:nx, :ny, :nlayers)
const XYLF = (:nx, :ny, :nlayers, :nfractions)

# (name, dims, unit, description)
const STATE_TABLE = Tuple{Symbol,Tuple,String,String}[(:x, XY, "m",
                                                       "Real-world x-coordinate of grid cell centres"),
                                                      (:y, XY, "m",
                                                       "Real-world y-coordinate of grid cell centres"),
                                                      (:ds, XY, "m", "Grid cell size in x"),
                                                      (:dn, XY, "m", "Grid cell size in y"),
                                                      (:dsdn, XY, "m^2", "Grid cell area"),
                                                      (:dsdni, XY, "m^-2",
                                                       "Inverse of grid cell area"),
                                                      (:zb, XY, "m",
                                                       "Bed level above reference"),
                                                      (:zs, XY, "m",
                                                       "Water level above reference (zb where zb > zs)"),
                                                      (:zne, XY, "m",
                                                       "Height of the non-erodible layer"),
                                                      (:zb0, XY, "m", "Initial bed level"),
                                                      (:zdry, XY, "m",
                                                       "Dry layer (unused)"),
                                                      (:dzdry, XY, "m",
                                                       "Dry layer change (unused)"),
                                                      (:dzb, XY, "m",
                                                       "Bed level change per time step"),
                                                      (:dzbyear, XY, "m/year",
                                                       "Bed level change per year"),
                                                      (:dzbavg, XY, "m/year",
                                                       "Bed level change averaged over avg_time"),
                                                      (:S, XY, "-", "Level of saturation"),
                                                      (:moist, XY, "-",
                                                       "Volumetric moisture content"),
                                                      (:moist_swr, XY, "-",
                                                       "Moisture content from soil water retention"),
                                                      (:h_delta, XY, "-",
                                                       "Suction at reversal between wetting and drying"),
                                                      (:gw, XY, "m", "Groundwater level"),
                                                      (:gw_prev, XY, "m",
                                                       "Groundwater level in the previous step"),
                                                      (:wetting, XY, "-",
                                                       "Wetting (1) or drying (0) of the soil profile"),
                                                      (:scan_w, XY, "-",
                                                       "Moisture on the wetting scanning curve (1/0)"),
                                                      (:scan_d, XY, "-",
                                                       "Moisture on the drying scanning curve (1/0)"),
                                                      (:scan_w_moist, XY, "-",
                                                       "Moisture on the wetting scanning curve"),
                                                      (:scan_d_moist, XY, "-",
                                                       "Moisture on the drying scanning curve"),
                                                      (:w_h, XY, "-",
                                                       "Moisture on the main wetting curve"),
                                                      (:d_h, XY, "-",
                                                       "Moisture on the main drying curve"),
                                                      (:w_hdelta, XY, "-",
                                                       "Moisture on the main wetting curve at h_delta"),
                                                      (:d_hdelta, XY, "-",
                                                       "Moisture on the main drying curve at h_delta"),
                                                      (:ustar, XY, "m/s", "Shear velocity"),
                                                      (:ustars, XY, "m/s",
                                                       "Shear velocity in x"),
                                                      (:ustarn, XY, "m/s",
                                                       "Shear velocity in y"),
                                                      (:ustar0, XY, "m/s",
                                                       "Shear velocity over a flat bed"),
                                                      (:ustars0, XY, "m/s",
                                                       "Shear velocity in x over a flat bed"),
                                                      (:ustarn0, XY, "m/s",
                                                       "Shear velocity in y over a flat bed"),
                                                      (:zsep, XY, "m",
                                                       "Level of the separation bubble"),
                                                      (:hsep, XY, "m",
                                                       "Height of the separation bubble above the bed"),
                                                      (:theta_dyn, XY, "deg",
                                                       "Dynamic angle of repose"),
                                                      (:rhoveg, XY, "-",
                                                       "Vegetation cover"),
                                                      (:drhoveg, XY, "-",
                                                       "Change in vegetation cover"),
                                                      (:hveg, XY, "m", "Vegetation height"),
                                                      (:dhveg, XY, "m/year",
                                                       "Vegetation height change"),
                                                      (:dzbveg, XY, "m/year",
                                                       "Bed level change for vegetation growth"),
                                                      (:germinate, XY, "-",
                                                       "Newly vegetated by germination (1/0)"),
                                                      (:lateral, XY, "-",
                                                       "Newly vegetated by lateral expansion (1/0)"),
                                                      (:vegetated, XY, "-",
                                                       "Vegetated (1/0)"),
                                                      (:vegfac, XY, "-",
                                                       "Shear reduction factor of vegetation (Raupach 1993)"),
                                                      (:fence_height, XY, "m",
                                                       "Fence height"),
                                                      (:fence_base, XY, "m",
                                                       "Bed level at the fences when placed"),
                                                      (:fence_top, XY, "m",
                                                       "Level of the fence tops"),
                                                      (:fence_height_init, XY, "m",
                                                       "Initial fence height"),
                                                      (:R, XY, "m", "Wave runup"),
                                                      (:eta, XY, "m", "Wave setup"),
                                                      (:sigma_s, XY, "m", "Swash"),
                                                      (:TWL, XY, "m",
                                                       "Total water level (SWL + runup)"),
                                                      (:SWL, XY, "m", "Still water level"),
                                                      (:DSWL, XY, "m",
                                                       "Dynamic still water level (SWL + setup)"),
                                                      (:Rti, XY, "-",
                                                       "Sheltering factor of roughness elements"),
                                                      (:gradh, XY, "-",
                                                       "Bed slope magnitude from avalanching"),
                                                      (:uw, XY, "m/s", "Wind velocity"),
                                                      (:uws, XY, "m/s",
                                                       "Wind velocity in x"),
                                                      (:uwn, XY, "m/s",
                                                       "Wind velocity in y"),
                                                      (:tau, XY, "N/m^2",
                                                       "Wind shear stress"),
                                                      (:taus, XY, "N/m^2",
                                                       "Wind shear stress in x"),
                                                      (:taun, XY, "N/m^2",
                                                       "Wind shear stress in y"),
                                                      (:tau0, XY, "N/m^2",
                                                       "Wind shear stress over a flat bed"),
                                                      (:taus0, XY, "N/m^2",
                                                       "Wind shear stress in x over a flat bed"),
                                                      (:taun0, XY, "N/m^2",
                                                       "Wind shear stress in y over a flat bed"),
                                                      (:taus_u, XY, "N/m^2",
                                                       "Saved direction of shear stress in x (unused)"),
                                                      (:taun_u, XY, "N/m^2",
                                                       "Saved direction of shear stress in y (unused)"),
                                                      (:dtaus, XY, "-",
                                                       "Shear perturbation in x (unused)"),
                                                      (:dtaun, XY, "-",
                                                       "Shear perturbation in y (unused)"),
                                                      (:udir, XY, "deg", "Wind direction"),
                                                      (:Hs, XY, "m", "Wave height"),
                                                      (:Hsmix, XY, "m",
                                                       "Wave height for mixing (including setup)"),
                                                      (:Tp, XY, "s", "Wave period"),
                                                      (:Cu, XYF, "kg/m^2",
                                                       "Equilibrium sediment concentration"),
                                                      (:Cuf, XYF, "kg/m^2",
                                                       "Equilibrium concentration with the fluid threshold"),
                                                      (:Cu0, XYF, "kg/m^2",
                                                       "Equilibrium concentration over a flat bed"),
                                                      (:Ct, XYF, "kg/m^2",
                                                       "Sediment concentration"),
                                                      (:q, XYF, "kg/m/s", "Sediment flux"),
                                                      (:qs, XYF, "kg/m/s",
                                                       "Sediment flux in x"),
                                                      (:qn, XYF, "kg/m/s",
                                                       "Sediment flux in y"),
                                                      (:pickup, XYF, "kg/m^2",
                                                       "Sediment entrainment"),
                                                      (:w, XYF, "-",
                                                       "Weights of sediment fractions"),
                                                      (:w_init, XYF, "-",
                                                       "Initial guess for w"),
                                                      (:w_air, XYF, "-",
                                                       "Weights from the grain size distribution in the air"),
                                                      (:w_bed, XYF, "-",
                                                       "Weights from the grain size distribution in the bed"),
                                                      (:uth, XYF, "m/s",
                                                       "Shear velocity threshold"),
                                                      (:uthf, XYF, "m/s",
                                                       "Fluid shear velocity threshold"),
                                                      (:uth0, XYF, "m/s",
                                                       "Shear velocity threshold from grain size only"),
                                                      (:u, XYF, "m/s",
                                                       "Saltation velocity"),
                                                      (:us, XYF, "m/s",
                                                       "Saltation velocity in x"),
                                                      (:un, XYF, "m/s",
                                                       "Saltation velocity in y"),
                                                      (:usST, XYF, "m/s",
                                                       "Saltation velocity in x for SedTRAILS"),
                                                      (:unST, XYF, "m/s",
                                                       "Saltation velocity in y for SedTRAILS"),
                                                      (:u0, XYF, "m/s",
                                                       "Saltation velocity over a flat bed"),
                                                      (:masstop, XYF, "kg/m^2",
                                                       "Sediment mass in the top layer"),
                                                      (:thlyr, XYL, "m",
                                                       "Bed layer thickness"),
                                                      (:salt, XYL, "-", "Salt content"),
                                                      (:mass, XYLF, "kg/m^2",
                                                       "Sediment mass in the bed")]

const STATE_NAMES = Tuple(first.(STATE_TABLE))
const STATE_INDEX = Dict(row[1] => row for row in STATE_TABLE)

let fields = [:($(row[1])::Array{T,$(length(row[2]))}) for row in STATE_TABLE]
    @eval begin
        """
            AeolisState{T}

        The spatial state of an AeoLiS model: one array per variable, of size `(nx+1, ny+1)`,
        `(nx+1, ny+1, nfractions)`, `(nx+1, ny+1, nlayers)` or
        `(nx+1, ny+1, nlayers, nfractions)`. Fields are read as `s.zb` or `s[:zb]`;
        [`state_info`](@ref) lists the variables with dimensions, units and descriptions.
        """
        mutable struct AeolisState{T<:AbstractFloat}
            $(fields...)
        end
    end
end

"""
    state_dims(name, nx, ny, nlayers, nfractions)

The size of the state variable `name` on a grid of `nx` by `ny` cells (`nx+1` by `ny+1`
points).
"""
function state_dims(name::Symbol, nx::Integer, ny::Integer, nl::Integer, nf::Integer)
    haskey(STATE_INDEX, name) || throw(ArgumentError("unknown state variable :$name"))
    d = Dict(:nx => nx + 1, :ny => ny + 1, :nlayers => nl, :nfractions => nf)
    return Tuple(d[k] for k in STATE_INDEX[name][2])
end

"""
    AeolisState{T}(nx, ny, nlayers, nfractions)

A state with every variable set to zero.
"""
function AeolisState{T}(nx::Integer, ny::Integer, nl::Integer, nf::Integer) where {T}
    arrays = [zeros(T, state_dims(r[1], nx, ny, nl, nf)) for r in STATE_TABLE]
    return AeolisState{T}(arrays...)
end

Base.getindex(s::AeolisState, name::Symbol) = getfield(s, name)
function Base.setindex!(s::AeolisState, v, name::Symbol)
    getfield(s, name) .= v
    return s
end
Base.keys(::AeolisState) = STATE_NAMES
Base.haskey(::AeolisState, name::Symbol) = haskey(STATE_INDEX, name)

function Base.copy(s::AeolisState{T}) where {T}
    return AeolisState{T}((copy(getfield(s, n)) for n in STATE_NAMES)...)
end

"""
    copyto_state!(dst::AeolisState, src::AeolisState)

Copy every variable of `src` into `dst`.
"""
function copyto_state!(dst::AeolisState, src::AeolisState)
    for n in STATE_NAMES
        copyto!(getfield(dst, n), getfield(src, n))
    end
    return dst
end

function Base.show(io::IO, ::MIME"text/plain", s::AeolisState{T}) where {T}
    nx, ny = size(s.zb)
    return print(io, "AeolisState{", T, "} on ", nx, " x ", ny, " points, ",
                 size(s.mass, 3), " layers, ", size(s.mass, 4), " fractions")
end
Base.show(io::IO, s::AeolisState{T}) where {T} = print(io, "AeolisState{", T, "}(...)")

"""
    state_info([name])

A table of every state variable as `(name, dims, unit, description)` named tuples, or the
entry of one variable.
"""
function state_info()
    return [(; name=r[1], dims=r[2], unit=r[3], description=r[4]) for r in STATE_TABLE]
end
function state_info(name::Symbol)
    haskey(STATE_INDEX, name) || throw(ArgumentError("unknown state variable :$name"))
    r = STATE_INDEX[name]
    return (; name=r[1], dims=r[2], unit=r[3], description=r[4])
end
