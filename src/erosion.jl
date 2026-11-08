# Wave-driven dune erosion (AeoLiS `erosion.py`): the model of Palmsten and Holman (2012) for
# collision of the total water level with the dune toe.

# The depth of sand cut from a cell, weighted by the resistance of its rooted layer.
function _resisted_depth(zb, zbT, zr, cr, tau0)
    d = zb - zbT
    rooted = clamp(zr, zero(d), d)
    return d + rooted * cr / tau0
end

_normcdf(x, mu, sigma) = 0.5 * SpecialFunctions.erfc(-(x - mu) / (sigma * sqrt(2.0)))

"""
    dune_erosion!(s, m, t)

Erode the dune face along each transect where the total water level (runup times 1.26 plus
the tide) exceeds `dune_toe_elevation`: the eroded volume `4 Cs (TWL - z_toe)^2 Nc`, with
`Nc` the number of collisions in the step, is removed by cutting the profile back to a slope
of `0.54 beach_slope` from the toe.

With marram grass (`veg_model = marram`), roots resist the cut: removing sand from a cell's
rooted layer (`veg_zr`) costs `1 + veg_cr / tau0` times as much of the wave demand as bare
sand, where `tau0` is the strength of moist bare sand ([Erikson2007](@citet)). Sand below the
rooted layer erodes as bare sand. For a uniform reinforcement deeper than the cut, this is the
erosion coefficient `Cs tau0 / (tau0 + c_r)`. Vegetated dunes that erode faster through
run-up and infiltration, as in the flume of [Feagin2023](@citet), are not represented.
"""
function dune_erosion!(s::AeolisState, m::AeolisModel, t::Real)
    p = m.p
    (m.tide !== nothing && m.wave_period !== nothing) ||
        throw(ConfigurationError("process_dune_erosion needs tide_file and wave_file"))
    zToe = p.dune_toe_elevation
    beta = p.beach_slope
    dt = m.dt_opt
    Kd = 1.26
    roots = p.process_vegetation && p.veg_model == :marram ? m.veg.traits : nothing
    wl = m.tide(t)
    Tp = m.wave_period(t)
    for iy in axes(s.zb, 2)
        twl = s.R[1, iy] * Kd + wl
        twl > zToe || continue
        x = view(s.x, :, iy)
        zb = view(s.zb, :, iy)
        eta, sigma_s = s.eta[1, iy], s.sigma_s[1, iy]
        dx = abs(x[2] - x[1])
        Bt = beta * 0.54
        st0 = argmin(abs.(zb .- zToe))
        zbT = fill(NaN, length(x))
        zbT[st0:end] .= Bt .* (x[st0:end] .- x[st0]) .+ zToe
        for k in eachindex(zbT)
            zbT[k] > zb[k] && (zbT[k] = zb[k])
        end
        Vc = roots === nothing ? cumsum(dx .* (zb[st0:end] .- zbT[st0:end])) :
             cumsum(dx .* _resisted_depth.(zb[st0:end], zbT[st0:end],
                                           view(s.veg_zr, st0:length(x), iy),
                                           view(s.veg_cr, st0:length(x), iy), roots.tau0))
        Vc .-= Vc[1]
        pc = 1 - _normcdf(zToe, eta + wl, sigma_s)
        Nc = pc * dt / Tp
        dV = 4 * p.Cs * (twl - zToe)^2 * Nc
        d = dV < 0 ? 0 : argmin(abs.(Vc .- dV)) - 1
        st = st0 + d
        zb[st0:(st - 1)] .= zbT[st0:(st - 1)]
    end
    return s
end
