# Wave-driven dune erosion (AeoLiS `erosion.py`): the model of Palmsten and Holman (2012) for
# collision of the total water level with the dune toe.

_normcdf(x, mu, sigma) = 0.5 * SpecialFunctions.erfc(-(x - mu) / (sigma * sqrt(2.0)))

"""
    dune_erosion!(s, m, t)

Erode the dune face along each transect where the total water level (runup times 1.26 plus
the tide) exceeds `dune_toe_elevation`: the eroded volume `4 Cs (TWL - z_toe)^2 Nc`, with
`Nc` the number of collisions in the step, is removed by cutting the profile back to a slope
of `0.54 beach_slope` from the toe.
"""
function dune_erosion!(s::AeolisState, m::AeolisModel, t::Real)
    p = m.p
    (m.tide !== nothing && m.wave_period !== nothing) ||
        throw(ConfigurationError("process_dune_erosion needs tide_file and wave_file"))
    zToe = p.dune_toe_elevation
    beta = p.beach_slope
    dt = m.dt_opt
    Kd = 1.26
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
        Vc = cumsum(dx .* (zb[st0:end] .- zbT[st0:end]))
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
