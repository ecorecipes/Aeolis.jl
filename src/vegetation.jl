# Vegetation (AeoLiS `vegetation.py`): cover and height, the reduction of shear by plants
# (Raupach 1993, or Okin 2008 behind each plant), germination, lateral expansion, growth and
# burial (the original AeoLiS growth or Duran and Moore 2013).

"""
    initialize_vegetation!(s, m)

Vegetation cover from `veg_file` (negative values set to zero), height
`hveg_max * sqrt(cover)`, and the vegetated cells. With `veg_model = marram`, marram grass is
set up from the cover ([`initialize_marram!`](@ref)).
"""
function initialize_vegetation!(s::AeolisState, m::AeolisModel)
    p = m.p
    if p.veg_file !== nothing
        s.rhoveg .= _as_grid(eltype(s.zb), p.veg_file, m.nx + 1, m.ny + 1, :veg_file)
        isnan(s.rhoveg[1, 1]) && (s.rhoveg .= 0)
    end
    s.rhoveg[s.rhoveg .< 0] .= 0
    s.hveg .= p.hveg_max .* sqrt.(s.rhoveg)
    s.vegetated .= s.rhoveg .> 0
    s.germinate .= 0
    s.lateral .= 0
    if p.vegshear_type == :okin
        m.okin = RotationGrid(p.dx, p.dy, 100.0)
    end
    p.process_vegetation && p.veg_model == :marram && initialize_marram!(s, m)
    return s
end

"""
    vegetation_shear!(s, m)

Reduce the shear velocity by vegetation: Raupach's factor `1 / sqrt(1 + gamma_vegshear *
cover)`, optionally smoothed by a Gaussian filter of width `veg_sigma`, or Okin's reduction
in the lee of each plant on a wind-aligned grid.
"""
function vegetation_shear!(s::AeolisState, m::AeolisModel)
    p = m.p
    if p.vegshear_type == :okin
        ustar = okin_shear!(m.okin, s, m, p.okin_c1_veg, p.okin_initialred_veg)
        s.ustar .= ustar
        @. s.ustars = -s.ustar * sin((-p.alfa + s.udir) / 180.0 * pi)
        @. s.ustarn = -s.ustar * cos((-p.alfa + s.udir) / 180.0 * pi)
    else
        ets = zeros(size(s.zb))
        etn = zeros(size(s.zb))
        for k in eachindex(s.ustar)
            if s.ustar[k] != 0
                ets[k] = s.ustars[k] / s.ustar[k]
                etn[k] = s.ustarn[k] / s.ustar[k]
            end
        end
        vegfac = 1.0 ./ sqrt.(1.0 .+ p.gamma_vegshear .* s.rhoveg)
        s.vegfac .= gaussian_filter(vegfac, p.veg_sigma)
        s.ustar .*= s.vegfac
        s.ustars .= s.ustar .* ets
        s.ustarn .= s.ustar .* etn
    end
    velocity_stress!(s, m)
    return s
end

"""
    germinate!(s, m)

Germination where the bed is not eroding, and lateral expansion from vegetated neighbours,
with yearly probabilities `germinate` and `lateral` scaled to the time step. Uses the
model's random number generator.
"""
function germinate!(s::AeolisState, m::AeolisModel)
    p = m.p
    dtacc = m.dt_opt * p.accfac
    n = SECONDS_PER_YEAR / dtacc
    pg = 1 - (1 - p.germinate)^(1.0 / n)
    rnd = rand(m.rng, size(s.zb)...)
    for k in eachindex(s.germinate)
        g = (s.dzbveg[k] >= 0) && (rnd[k] <= pg)
        s.germinate[k] = g
        s.vegetated[k] = (g || s.vegetated[k] > 0)
    end
    dx = m.ny > 1 ? s.ds[3, 3] : p.dx
    pl = 1 - (1 - p.lateral)^(1.0 / n)
    plc = 1 - (1 - pl)^(1.0 / dx)
    rho, ds, dn = s.rhoveg, s.ds, s.dn
    nx, ny = size(rho)
    for j in 1:ny, i in 1:nx
        c = 0
        i > 1 && max((rho[i - 1, j] - rho[i, j]) / ds[i, j], 0.0) > 0 && (c += 1)
        i < nx && max((rho[i + 1, j] - rho[i, j]) / ds[i, j], 0.0) > 0 && (c += 1)
        j > 1 && max((rho[i, j - 1] - rho[i, j]) / dn[i, j], 0.0) > 0 && (c += 1)
        j < ny && max((rho[i, j + 1] - rho[i, j]) / dn[i, j], 0.0) > 0 && (c += 1)
        s.drhoveg[i, j] = c
    end
    for k in eachindex(s.vegetated)
        s.vegetated[k] = (rnd[k] <= plc * s.drhoveg[k]) || s.vegetated[k] > 0
    end
    return s
end

"""
    grow_vegetation!(s, m)

Vertical growth towards `hveg_max` at rate `V_ver`, reduced by burial or erosion away from
`dzb_opt`; plants are removed where flooded or below `veg_min_elevation`.
"""
function grow_vegetation!(s::AeolisState, m::AeolisModel)
    p = m.p
    dtacc = m.dt_opt * p.accfac
    ix = (s.vegetated .> 0) .& (p.V_ver > 0)
    vgrow = zeros(size(s.hveg))
    vgrow[ix] .= p.V_ver .* (1 .- s.hveg[ix] ./ p.hveg_max)
    apply_mask!(vgrow, m.masks[:vver_mask])
    burial = abs.(s.dzbveg .- p.dzb_opt) .* p.veg_gamma
    s.dhveg[ix] .= vgrow[ix] .- burial[ix]
    if p.veggrowth_type == :orig
        s.hveg .+= s.dhveg .* dtacc ./ SECONDS_PER_YEAR
        s.hveg .= max.(min.(s.hveg, p.hveg_max), 0.0)
        s.rhoveg .= (s.hveg ./ p.hveg_max) .^ 2
    else
        t_veg = p.t_veg / 365
        s.rhoveg[s.rhoveg .> p.rhoveg_max] .= p.rhoveg_max
        ixzero = s.rhoveg .<= 0
        if p.V_ver > 0
            # AeoLiS stores drhoveg as the integer neighbour count from `germinate`, so the
            # growth rate assigned to it is truncated towards zero
            s.drhoveg[ix] .= trunc.((p.rhoveg_max .- s.rhoveg[ix]) ./ t_veg .-
                                    (p.v_gam / p.hveg_max) .*
                                    abs.(s.dzbveg[ix] .- p.dzb_opt) .*
                                    p.veg_gamma)
        else
            s.drhoveg[ix] .= 0
        end
        s.rhoveg .+= s.drhoveg .* (p.dt * p.accfac) ./ SECONDS_PER_YEAR
        s.rhoveg[s.rhoveg .< 0] .= 0
        s.rhoveg[ixzero] .= 0
        s.hveg .= p.hveg_max .* sqrt.(s.rhoveg)
    end
    dhmax = p.V_ver * dtacc / SECONDS_PER_YEAR
    drmax = (dhmax / p.hveg_max)^2
    s.vegetated .*= s.rhoveg .>= 0.5 * drmax
    if p.process_tide
        fl = s.zb .< s.TWL
        s.rhoveg[fl] .= 0
        s.hveg[fl] .= 0
        s.vegetated[fl] .= 0
    end
    low = s.zb .< p.veg_min_elevation
    s.rhoveg[low] .= 0
    s.hveg[low] .= 0
    s.vegetated[low] .= 0
    return s
end

# Okin's shear reduction along each wind-aligned row (AeoLiS `compute_okin_shear`).
function _okin_rows(x::AbstractMatrix, ustar::AbstractMatrix, hveg::AbstractMatrix, c1,
                    intercept)
    intercept = min(intercept, 1.0)
    nx, ny = size(x)
    out = copy(ustar)
    for j in 1:ny
        red_all = zeros(nx)
        for g in 1:nx
            h = hveg[g, j]
            h > 0 || continue
            for k in 1:nx
                xrel = x[k, j] - x[g, j]
                mult = xrel >= 0 ? intercept + (1 - intercept) * (1 - exp(-xrel * c1 / h)) :
                       1.0
                red = 1 - mult
                red < 0.00001 && (red = 0.0)
                red > 1 && (red = 1.0)
                xrel < 0 && (red = 0.0)
                red_all[k] += red
            end
        end
        for k in 1:nx
            r = min(red_all[k], 1.0)
            out[k, j] = max(ustar[k, j] * (1 - r), 0.01)
        end
    end
    return out
end
