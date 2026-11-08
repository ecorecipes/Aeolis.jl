# Marram grass on the grid (`veg_model = marram`): every cell is a `MarramCell` held in the
# `veg_*` state variables. Each time step, the canopy sets the effective cover and height
# that the existing shear reduction uses; every `veg_interval` of accelerated time, every cell
# is advanced by `marram_step` and bare cells are colonised from their neighbours.

"""
    initialize_marram!(s, m)

Set up marram grass on the grid:
- the traits from `veg_traits_file` (or the defaults), and the rainfall and potential
  evapotranspiration of `veg_water_file`;
- the initial shoot density from `veg_file`, read as a fraction of `N_max`;
- the canopy at `h_new`, roots at their target without burial, and a rooted layer of
  `zr_new`;
- root-zone moisture at field capacity, no soil-biota pressure;
- a colonisation threshold for every cell, from the model's random number generator.
"""
function initialize_marram!(s::AeolisState, m::AeolisModel)
    p = m.p
    tr = marram_traits(p)
    w = p.veg_water_file
    rain = w === nothing ? nothing : CircularSeries(w[:, 1], w[:, 2])
    pet = w === nothing ? nothing : CircularSeries(w[:, 1], w[:, 3])
    m.veg = (; traits=tr, rain, pet)
    frac = clamp.(s.rhoveg, 0, 1)
    @. s.veg_N = ifelse(frac * tr.N_max >= tr.N_min, frac * tr.N_max, 0.0)
    for k in eachindex(s.veg_N)
        s.veg_Br[k] = root_target(s.veg_N[k], 0.0, tr)
    end
    s.veg_zr .= ifelse.(s.veg_N .> 0, tr.zr_new, 0.0)
    s.veg_s .= tr.s_fc
    s.veg_P .= 0
    s.veg_Dsub .= 0
    s.veg_tdry .= 0
    s.veg_zdot .= 0
    s.veg_zbref .= s.zb
    s.veg_tacc .= 0
    s.veg_wet .= 0
    s.veg_H .= 0
    _draw_thresholds!(s, m, trues(size(s.veg_E)))
    s.hveg .= ifelse.(s.veg_N .> 0, tr.h_new, 0.0)
    s.vegetated .= s.veg_N .> 0
    marram_cover!(s, m)
    return s
end

_is_transect(m::AeolisModel) = ndims(m.p.xgrid_file) == 1

# New colonisation thresholds where `mask`: exponential draws, or 1 (the expected value) when
# `veg_random` is off. A transect draws one value per column, so its three rows stay equal.
function _draw_thresholds!(s::AeolisState, m::AeolisModel, mask)
    E = s.veg_E
    if !m.p.veg_random
        E[mask] .= 1
    elseif _is_transect(m)
        for i in axes(E, 1)
            any(@view mask[i, :]) || continue
            e = Random.randexp(m.rng)
            for j in axes(E, 2)
                mask[i, j] && (E[i, j] = e)
            end
        end
    else
        for k in eachindex(E)
            mask[k] && (E[k] = Random.randexp(m.rng))
        end
    end
    return s
end

"""
    marram_cover!(s, m)

The vegetation cover and height seen by the shear reduction: `hveg` is the canopy height and
`rhoveg` the [`effective_cover`](@ref) of the frontal area index, which depends on the
current wind speed through the reconfiguration of the shoots ([`frontal_area`](@ref)).
"""
function marram_cover!(s::AeolisState, m::AeolisModel)
    tr = m.veg.traits
    g = m.p.gamma_vegshear
    for k in eachindex(s.veg_N)
        lam = frontal_area(s.veg_N[k], s.hveg[k], s.uw[k], tr)
        s.veg_lam[k] = lam
        s.rhoveg[k] = g > 0 ? effective_cover(lam, tr, g) : 0.0
    end
    return s
end

"""
    marram_update!(s, m, t)

Advance marram grass after the time step that started at `t`. The time step, accelerated by
`accfac`, is added to the time since the last update, and to the time under water of cells
whose water level is above the bed. Once `veg_interval` has passed, every cell is advanced by
[`marram_step`](@ref), with:
- the mean rate of bed-level change since the last update;
- the day of year at the middle of the interval, from `refdate` and accelerated time;
- rainfall and potential evapotranspiration from `veg_water_file`, at the same time;
- the fraction of the interval under water;
- the current wind speed.

With groundwater, a root zone whose base is within 1 m of the water table is kept at field
capacity, the capillary supply of [Homberger2024](@citet). Below `veg_min_elevation`, plants
are removed, as in AeoLiS. Bare cells are then colonised ([`marram_colonise!`](@ref)).
"""
function marram_update!(s::AeolisState, m::AeolisModel, t::Real)
    p = m.p
    tr = m.veg.traits
    dtacc = m.dt_opt * p.accfac
    s.veg_tacc .+= dtacc
    if p.process_tide || p.process_wave || p.process_runup
        for k in eachindex(s.veg_wet)
            s.zs[k] > s.zb[k] && (s.veg_wet[k] += dtacc)
        end
    end
    interval = s.veg_tacc[1]
    interval < p.veg_interval && return s
    days = interval / 86400
    tmid = (t + m.dt_opt) * p.accfac - interval / 2
    doy = day_of_year(model_datetime(p, tmid))
    rain = m.veg.rain === nothing ? nothing : m.veg.rain(tmid)
    pet = m.veg.pet === nothing ? 0.0 : m.veg.pet(tmid)
    for k in eachindex(s.veg_N)
        c = MarramCell(s.veg_N[k], s.hveg[k], s.veg_Br[k], s.veg_zr[k], s.veg_s[k],
                       s.veg_P[k], s.veg_Dsub[k], s.veg_tdry[k], s.veg_zdot[k])
        dzb = (s.zb[k] - s.veg_zbref[k]) / days
        wet = s.veg_wet[k] / interval
        f = rain === nothing ? (; doy, dzb, u=s.uw[k], submerged=wet) :
            (; doy, dzb, u=s.uw[k], submerged=wet, rain, pet)
        c = marram_step(c, f, tr, days)
        sw = c.s
        if p.process_groundwater && s.zb[k] - s.gw[k] < tr.Z_r + 1
            sw = max(sw, tr.s_fc)
        end
        N, h = c.N, c.h
        if s.zb[k] < p.veg_min_elevation
            N, h = 0.0, 0.0
        end
        s.veg_N[k] = N
        s.hveg[k] = h
        s.veg_Br[k] = c.Br
        s.veg_zr[k] = c.zr
        s.veg_s[k] = sw
        s.veg_P[k] = c.P
        s.veg_Dsub[k] = c.Dsub
        s.veg_tdry[k] = c.tdry
        s.veg_zdot[k] = c.zdot
        r = marram_rates(MarramCell(N, h, c.Br, c.zr, sw, c.P, c.Dsub, c.tdry, c.zdot), f,
                         tr)
        s.veg_Gb[k] = r.Gb
        s.veg_fW[k] = r.fW
        s.veg_V[k] = r.V
        s.veg_growth[k] = r.growth
        s.veg_mort[k] = r.mortality
        s.veg_cr[k] = r.cr
    end
    marram_colonise!(s, m, season(doy, tr), days)
    s.vegetated .= s.veg_N .> 0
    s.veg_zbref .= s.zb
    s.veg_tacc .= 0
    s.veg_wet .= 0
    marram_cover!(s, m)
    return s
end

"""
    marram_colonise!(s, m, phi, days)

Colonisation of bare cells over `days` days with season weight `phi`. A bare cell
accumulates the hazard `phi (v_front / L_eff) sum(w o_nb) / dx` per day, plus `establishment`
per year:
- `L_eff = L_s - season_ramp` is the number of full growing days in a year;
- `o_nb = min(N_nb / (N_new N_max), 1)` is the occupancy of a neighbour. Rhizomes extend from
  any established clone, however sparse, so a front advances at about `v_front` whatever the
  grid spacing;
- the sum runs over its neighbours along each axis (only along x on a transect), with
  weight `1 - f_jump`, and over the cells `jump_distance` (in m) away, with weight `f_jump`.

When the accumulated hazard reaches the cell's threshold (exponential, or 1 when `veg_random`
is off), the cell is colonised with `N_new N_max` shoots, a canopy of `h_new` and a rooted
layer of `zr_new`. The front speed `v_front` is that of patches on Ameland
([Keijsers2015](@citet)), and the jumps those of [Bonte2021](@citet).

Cells below `veg_min_elevation`, and cells that were under water for more than half the
interval, are not colonised.
"""
function marram_colonise!(s::AeolisState, m::AeolisModel, phi::Real, days::Real)
    p = m.p
    tr = m.veg.traits
    N = s.veg_N
    nx, ny = size(N)
    transect = _is_transect(m)
    a = tr.v_front / (season_length(tr) - tr.season_ramp)
    occ(n) = min(n / (tr.N_new * tr.N_max), 1.0)
    estab = tr.establishment / DAYS_PER_YEAR
    rows = transect ? (2:2) : (1:ny)
    newly = falses(size(N))
    for j in rows, i in 1:nx
        if N[i, j] > 0
            s.veg_H[i, j] = 0
            continue
        end
        (s.zb[i, j] < p.veg_min_elevation || s.veg_wet[i, j] > 0.5 * s.veg_tacc[i, j]) &&
            continue
        cx = 0.0
        di_jump = max(round(Int, tr.jump_distance / s.dn[i, j]), 1)
        for (di, w) in ((1, 1 - tr.f_jump), (di_jump, tr.f_jump))
            i > di && (cx += w * occ(N[i - di, j]))
            i + di <= nx && (cx += w * occ(N[i + di, j]))
        end
        cy = 0.0
        if !transect
            dj_jump = max(round(Int, tr.jump_distance / s.ds[i, j]), 1)
            for (dj, w) in ((1, 1 - tr.f_jump), (dj_jump, tr.f_jump))
                j > dj && (cy += w * occ(N[i, j - dj]))
                j + dj <= ny && (cy += w * occ(N[i, j + dj]))
            end
        end
        hz = phi * (a * (cx / s.dn[i, j] + cy / s.ds[i, j]) + estab)
        s.veg_H[i, j] += hz * days
        s.veg_H[i, j] >= s.veg_E[i, j] && (newly[i, j] = true)
    end
    if transect
        for j in (1, 3)
            newly[:, j] .= view(newly, :, 2)
            s.veg_H[:, j] .= view(s.veg_H, :, 2)
        end
    end
    for k in eachindex(N)
        newly[k] || continue
        N[k] = tr.N_new * tr.N_max
        s.hveg[k] = tr.h_new
        s.veg_zr[k] = tr.zr_new
        s.veg_Br[k] = max(s.veg_Br[k], root_target(N[k], 0.0, tr))
        s.veg_H[k] = 0
    end
    _draw_thresholds!(s, m, newly)
    return s
end
