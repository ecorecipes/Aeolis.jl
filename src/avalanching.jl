# Avalanching (AeoLiS `avalanching.py`): sand moves down slopes steeper than the dynamic angle
# of repose until no slope exceeds it.

"""
    angle_of_repose!(s, m)

Set the dynamic angle of repose of every cell to `theta_dyn`.
"""
function angle_of_repose!(s::AeolisState, m::AeolisModel)
    s.theta_dyn .= m.p.theta_dyn
    return s
end

"""
    avalanche!(s, m)

Redistribute sand down slopes steeper than the angle of repose (an explicit iteration with a
relaxation factor 0.1, at most `max_iter_ava` sweeps; the downslope gradient of a cell is the
steeper one-sided difference, zero on the grid edges and at local minima; fluxes wrap
around the grid). The water level is reset to at least the bed; like AeoLiS, the still water
level then takes the same values.
"""
function avalanche!(s::AeolisState, m::AeolisModel)
    p = m.p
    p.process_avalanche || return s
    tan_dyn = tan(deg2rad(p.theta_dyn))
    zb = s.zb
    nx, ny = size(zb)
    E = 0.1
    gx = zeros(nx, ny)
    gy = zeros(nx, ny)
    fx = zeros(nx, ny)
    fy = zeros(nx, ny)
    gh = zeros(nx, ny)
    for _ in 1:(p.max_iter_ava)
        fill!(gx, 0)
        fill!(gy, 0)
        fill!(fx, 0)
        fill!(fy, 0)
        fill!(gh, 0)
        maxg = 0.0
        @inbounds for j in 1:ny, i in 1:nx
            s.zne[i, j] >= zb[i, j] && continue
            c = zb[i, j]
            g0 = 0.0
            if !(i == 1 || i == nx)
                r, l = zb[i + 1, j], zb[i - 1, j]
                if !((r > c) && (l > c))
                    g0 = r > l ? l - c : c - r
                end
                gx[i, j] = g0 / s.ds[i, j]
            end
            g1 = 0.0
            if !(j == 1 || j == ny)
                d, u = zb[i, j + 1], zb[i, j - 1]
                if !((d > c) && (u > c))
                    g1 = d > u ? u - c : c - d
                end
                gy[i, j] = g1 / s.dn[i, j]
            end
            g = sqrt(gx[i, j]^2 + gy[i, j]^2)
            gh[i, j] = g
            g > maxg && (maxg = g)
        end
        maxg < tan_dyn && break
        @inbounds for j in 1:ny, i in 1:nx
            if gh[i, j] > tan_dyn
                sd = tanh(gh[i, j]) - tanh(0.9 * tan_dyn)
                fx[i, j] = sd * gx[i, j]
                fy[i, j] = sd * gy[i, j]
            end
        end
        @inbounds for j in 1:ny, i in 1:nx
            qout = max(fx[i, j], 0.0) + max(-fx[i, j], 0.0) + max(fy[i, j], 0.0) +
                   max(-fy[i, j], 0.0)
            im = i == 1 ? nx : i - 1
            ip = i == nx ? 1 : i + 1
            jm = j == 1 ? ny : j - 1
            jp = j == ny ? 1 : j + 1
            qin = max(fx[im, j], 0.0) + max(-fx[ip, j], 0.0) + max(fy[i, jm], 0.0) +
                  max(-fy[i, jp], 0.0)
            zb[i, j] += E * (qin - qout)
        end
    end
    s.gradh .= gh
    for k in eachindex(s.SWL)
        zb[k] > s.SWL[k] && (s.SWL[k] = zb[k])
    end
    s.zs .= s.SWL
    return s
end
