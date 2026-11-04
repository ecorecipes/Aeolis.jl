# Equilibrium sediment transport (AeoLiS `transport.py`): grain velocities (wind speed, Duran
# et al. 2011, or constant), the saturated concentration by several formulations, and the
# weights of the fractions in the air and in the bed.

"""
    grainspeed_duran(s, m)

Grain velocities after Duran et al. (2011), with the effect of bed slopes; also the velocity
over a flat bed and the SedTRAILS velocities (zero where the wind is below the threshold).
Returns `(u0, us, un, u, usST, unST)`, each `(nx+1, ny+1, nfractions)`.
"""
function grainspeed_duran(s::AeolisState, m::AeolisModel)
    p = m.p
    nf = m.nf
    d = m.grain_size
    z, x, y = s.zb, s.x, s.y
    srho = p.rhog / p.rhoa
    kappa = p.kappa
    g, v = p.g, p.v
    sz = size(s.uth0)
    u0 = zeros(sz)
    us = zeros(sz)
    un = zeros(sz)
    u = zeros(sz)
    usST = zeros(sz)
    unST = zeros(sz)
    dzs = zeros(size(z))
    dzn = zeros(size(z))
    dzs[2:(end - 1), :] .= (z[3:end, :] .- z[1:(end - 2), :]) ./
                           (x[3:end, :] .- x[1:(end - 2), :])
    dzn[:, 2:(end - 1)] .= (z[:, 1:(end - 2)] .- z[:, 3:end]) ./
                           (y[:, 1:(end - 2)] .- y[:, 3:end])
    r = 1.0
    c = 14.0 / (1.0 + 1.4 * r)
    tv = (v / g^2)^(1 / 3)
    lv = (v^2 / (p.Aa^2 * g * (srho - 1)))^(1 / 3)
    z1 = 35.0 * lv
    A, B = 0.95, 5.12
    for q in 1:nf
        z0 = d[q] / 20.0
        alpha = 0.17 * d[q] / lv
        Sstar = d[q] / (4 * v) * sqrt(g * d[q] * (srho - 1.0))
        Cd = (4 / 3) * (A + sqrt(2 * alpha) * B / Sstar)^2
        uf = sqrt(4 / (3 * Cd) * (srho - 1.0) * g * d[q])
        for j in axes(z, 2), i in axes(z, 1)
            uth = s.uth0[i, j, q]
            zm = c * uth * tv
            ets, etn = 1.0, 0.0
            if s.uw[i, j] > 0
                ets = s.uws[i, j] / s.uw[i, j]
                etn = s.uwn[i, j] / s.uw[i, j]
            end
            if s.ustar[i, j] > 0
                ets = s.ustars[i, j] / s.ustar[i, j]
                etn = s.ustarn[i, j] / s.ustar[i, j]
            end
            dhs, dhn = dzs[i, j], dzn[i, j]
            Axs = ets + 2 * alpha * dhs
            Axn = etn + 2 * alpha * dhn
            Ax = hypot(Axs, Axn)
            ueff = (uth / kappa) * (log(z1 / z0) + (z1 / zm) * (s.ustar[i, j] / uth - 1))
            ueff0 = (uth / kappa) * (log(z1 / z0) + (z1 / zm) * (s.ustar0[i, j] / uth - 1))
            u0[i, j, q] = ueff0 - uf / sqrt(2 * alpha)
            usa = (ueff - uf / (sqrt(2.0 * alpha) * Ax)) * ets -
                  (sqrt(2 * alpha) * uf / Ax) * dhs
            una = (ueff - uf / (sqrt(2.0 * alpha) * Ax)) * etn -
                  (sqrt(2 * alpha) * uf / Ax) * dhn
            if p.method_grainspeed == :duran
                us[i, j, q], un[i, j, q] = usa, una
                u[i, j, q] = hypot(usa, una)
            else # :duran_full, Newton iterations on the grain momentum balance
                ui = _duran_full(complex(usa, una), complex(ueff * ets, ueff * etn), uf,
                                 alpha, complex(dhs, dhn))
                us[i, j, q], un[i, j, q] = real(ui), imag(ui)
                u[i, j, q] = abs(ui)
            end
            usST[i, j, q] = us[i, j, q]
            unST[i, j, q] = un[i, j, q]
            if s.uth[i, j, q] > s.ustar[i, j]
                usST[i, j, q] *= 0
                unST[i, j, q] *= 0
            end
        end
    end
    return u0, us, un, u, usST, unST
end

# scipy.optimize.newton with a secant method (no derivative), maxiter 20, tol 0.05
function _duran_full(u0::Complex, veff::Complex, uf, alpha, dh::Complex)
    f(ui) = (veff - ui) * abs(veff - ui) / uf^2 - ui / (2 * alpha * abs(ui)) - dh
    p0 = u0
    p1 = u0 * (1 + 1e-4) + (real(u0) >= 0 ? 1e-4 : -1e-4)
    q0, q1 = f(p0), f(p1)
    abs(q1) < abs(q0) && ((p0, p1, q0, q1) = (p1, p0, q1, q0))
    for _ in 1:20
        q1 == q0 && return (p1 + p0) / 2
        pn = p1 - q1 * (p1 - p0) / (q1 - q0)
        abs(pn - p1) < 0.05 && return pn
        p0, q0 = p1, q1
        p1, q1 = pn, f(pn)
    end
    return p1
end

"""
    grainspeed_constant!(s, m)

Grain velocity of 1 m/s along the shear velocity (or along the wind where the shear velocity
is zero).
"""
function grainspeed_constant!(s::AeolisState, m::AeolisModel)
    for j in axes(s.zb, 2), i in axes(s.zb, 1)
        us = un = 0.0
        if s.ustar[i, j] != 0
            us = s.ustars[i, j] / s.ustar[i, j]
            un = s.ustarn[i, j] / s.ustar[i, j]
        end
        if s.ustar[i, j] == 0 && s.uw[i, j] != 0
            us = s.uws[i, j] / s.uw[i, j]
            un = s.uwn[i, j] / s.uw[i, j]
        end
        s.us[i, j, :] .= us
        s.un[i, j, :] .= un
        s.u[i, j, :] .= hypot(us, un)
    end
    return s
end

"""
    equilibrium!(s, m)

Grain velocities and the saturated sediment concentration `Cu` (and `Cuf` with the fluid
threshold, `Cu0` over a flat bed) by `method_transport`: `:bagnold`, `:bagnold_gs`,
`:kawamura`, `:lettau`, `:dk`, `:sauermann` or `:vanrijn_strypsteen`, multiplied by `accfac`.
"""
function equilibrium!(s::AeolisState, m::AeolisModel)
    p = m.p
    if p.process_transport
        nf = m.nf
        if p.method_grainspeed in (:duran, :duran_full)
            u0, us, un, u, usST, unST = grainspeed_duran(s, m)
            s.u0 .= u0
            s.us .= us
            s.un .= un
            s.u .= u
            s.usST .= usST
            s.unST .= unST
        elseif p.method_grainspeed == :windspeed
            for q in 1:nf
                s.u0[:, :, q] .= s.uw
                s.us[:, :, q] .= s.uws
                s.un[:, :, q] .= s.uwn
                s.u[:, :, q] .= s.uw
            end
        else
            grainspeed_constant!(s, m)
        end
        rhoa, g = p.rhoa, p.g
        meth = p.method_transport
        s.Cu .= 0
        s.Cuf .= 0
        for q in 1:nf, j in axes(s.zb, 2), i in axes(s.zb, 1)
            ust = s.ustar[i, j]
            ust0 = s.ustar0[i, j]
            uu = s.u[i, j, q]
            (ust != 0 && uu != 0) || continue
            uth, uthf, uth0 = s.uth[i, j, q], s.uthf[i, j, q], s.uth0[i, j, q]
            if meth == :bagnold
                s.Cu[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust - uth)^3 / uu)
                s.Cuf[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust - uthf)^3 / uu)
                s.Cu0[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust0 - uth0)^3 / uu)
            elseif meth == :bagnold_gs
                f = p.Cb * sqrt(m.grain_size[q] / 0.000250) * rhoa / g * (ust - uth)^3 / uu
                s.Cu[i, j, q] = max(0.0, f)
                s.Cuf[i, j, q] = max(0.0, f)
                s.Cu0[i, j, q] = max(0.0, f)
            elseif meth == :kawamura
                s.Cu[i, j, q] = max(0.0, p.Ck * rhoa / g * (ust + uth)^2 * (ust - uth) / uu)
                s.Cuf[i, j, q] = max(0.0,
                                     p.Ck * rhoa / g * (ust + uthf)^2 * (ust - uthf) / uu)
            elseif meth == :lettau
                s.Cu[i, j, q] = max(0.0, p.Cl * rhoa / g * ust^2 * (ust - uth) / uu)
                s.Cuf[i, j, q] = max(0.0, p.Cl * rhoa / g * ust^2 * (ust - uthf) / uu)
            elseif meth == :dk
                s.Cu[i, j, q] = max(0.0,
                                    p.Cdk * rhoa / g * 0.8 * uth * (ust^2 - (0.8 * uth)^2) /
                                    uu)
                s.Cuf[i, j, q] = max(0.0,
                                     p.Cdk * rhoa / g * 0.8 * uthf *
                                     (ust^2 - (0.8 * uthf)^2) / uu)
                s.Cu0[i, j, q] = max(0.0,
                                     p.Cdk * rhoa / g * 0.8 * uth0 *
                                     (ust0^2 - (0.8 * uth0)^2) / uu)
            elseif meth == :sauermann
                a = 0.35
                s.Cu[i, j, q] = max(0.0, 2.0 * a * rhoa / g * (ust^2 - uth^2))
                s.Cuf[i, j, q] = max(0.0, 2.0 * a * rhoa / g * (ust^2 - uthf^2))
                s.Cu0[i, j, q] = max(0.0, 2.0 * a * rhoa / g * (ust0^2 - uth0^2))
            else # :vanrijn_strypsteen
                s.Cu[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust^3 - uth^3) / uu)
                s.Cuf[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust^3 - uth^3) / uu)
                s.Cu0[i, j, q] = max(0.0, p.Cb * rhoa / g * (ust0^3 - uth0^3) / uu)
            end
        end
    end
    s.Cu .*= p.accfac
    s.Cuf .*= p.accfac
    s.Cu0 .*= p.accfac
    return s
end

"""
    compute_weights(s, m)

Weights of the fractions: from the concentrations in the air relative to saturation
(`w_air`), from the top layer of the bed (`w_bed`), and their combination by the bed
interaction factor `bi`, normalised. Returns `(w, w_air, w_bed)`.
"""
function compute_weights(s::AeolisState, m::AeolisModel)
    bi = m.p.bi
    w_air = normalize_along(s.Ct, 3; ref=s.Cu)
    w_bed = normalize_along(s.mass[:, :, 1, :], 3)
    sa = sum(w_air; dims=3)
    w = (1.0 - bi) .* w_air .+ (1.0 .- min.(1.0, (1.0 - bi) .* sa)) .* w_bed
    return normalize_along(w, 3), w_air, w_bed
end

"""
    renormalize_weights!(w, ix)

Rescale the weights of fractions `ix` and above so that all weights sum to one; cells whose
weights are all zero get equal weights.
"""
function renormalize_weights!(w::AbstractArray, ix::Integer)
    nf = size(w, 3)
    f = sum(w[:, :, 1:(ix - 1)]; dims=3)
    rest = normalize_along(w[:, :, ix:nf], 3)
    w[:, :, ix:nf] .= rest .* (1.0 .- f)
    s = sum(w; dims=3)
    for k in CartesianIndices(w)
        sk = s[k[1], k[2], 1]
        w[k] = sk != 0 ? w[k] / sk : 1.0 / nf
    end
    return w
end
