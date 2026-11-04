# Wind (AeoLiS `wind.py`): the wind at the current time, the roughness length, the shear
# velocity of the undisturbed wind and its conversion to and from shear stress, and the
# dispatch to the shear perturbation methods.

"""
    initialize_wind!(s, m)

Prepare the shear computation for two-dimensional grids (transects are quasi-2-D, so they
count): the spectral `WindShear` for `method_shear = :fft`, a rotated grid for `:quasi2d` and
`:duna2d`.
"""
function initialize_wind!(s::AeolisState, m::AeolisModel)
    p = m.p
    z0 = roughness_length(s, m)
    m.z0_init = z0
    if p.process_shear && m.ny > 0
        if p.method_shear == :fft
            m.shear = WindShear(p.dx, p.dy, p.L, p.l, z0, p.buffer_width)
        else
            m.shear = RotationGrid(p.dx, p.dy, 100.0)
        end
    end
    return s
end

"""
    interpolate_wind!(s, m, t)

Wind speed and direction at time `t` (the value at the last time not after `t`, repeating the
series beyond its end; directions through their sine and cosine), its components on the grid,
and the shear velocity and stress of the undisturbed wind over a bed of roughness `z0`.
"""
function interpolate_wind!(s::AeolisState, m::AeolisModel, t::Real)
    p = m.p
    if p.process_wind && m.wind !== nothing
        w = m.wind
        k = findfirst(==(t), view(w, :, 1))
        if k !== nothing
            s.uw .= w[k, 2]
            s.udir .= w[k, 3]
        else
            s.uw .= m.wind_speed(t)
            s.udir .= atan(m.wind_sin(t), m.wind_cos(t)) * 180.0 / pi
        end
    end
    a = p.alfa
    @. s.uws = -s.uw * sin((-a + s.udir) / 180.0 * pi)
    @. s.uwn = -s.uw * cos((-a + s.udir) / 180.0 * pi)
    s.uw .= abs.(s.uw)
    z0 = roughness_length(s, m)
    @. s.ustars = s.uws * p.kappa / log(p.z / z0)
    @. s.ustarn = s.uwn * p.kappa / log(p.z / z0)
    @. s.ustar = hypot(s.ustars, s.ustarn)
    velocity_stress!(s, m)
    s.ustar0 .= s.ustar
    s.ustars0 .= s.ustars
    s.ustarn0 .= s.ustarn
    s.tau0 .= s.tau
    s.taus0 .= s.taus
    s.taun0 .= s.taun
    return s
end

"""
    roughness_length(s, m)

The roughness length `z0` for `method_roughness`: `:constant` (`z0 = k`, as AeoLiS for
backward compatibility), `:constant_nikuradse` (`k/30`), `:mean_grainsize_initial`,
`:mean_grainsize_adaptive`, `:median_grainsize_adaptive` (Sherman and Greenwood 1982) or
`:vanrijn_strypsteen` (van Rijn and Strypsteen 2019). A scalar or a grid.
"""
function roughness_length(s::AeolisState, m::AeolisModel)
    p = m.p
    meth = p.method_roughness
    if meth == :constant
        return p.k
    elseif meth == :constant_nikuradse
        return p.k / 30
    elseif meth == :mean_grainsize_initial
        return sum(p.grain_size .* vec(m.grain_dist[1, :])) / 30.0
    elseif meth == :mean_grainsize_adaptive
        return mean_grain_size(m.grain_size, s.mass) ./ 30.0
    elseif meth == :median_grainsize_adaptive
        return 2 .* grain_size_percentile(m.grain_size, s.mass, 50) ./ 30.0
    else # :vanrijn_strypsteen
        if length(p.grain_dist) == 1
            d50 = p.grain_size[1]
            d90 = 2 * d50
        else
            d50 = grain_size_percentile(m.grain_size, s.mass, 50)
            d90 = grain_size_percentile(m.grain_size, s.mass, 90)
        end
        ustar_grain = @. p.kappa * (s.uw / log(30 * p.z / d90))
        ustar_th_B = @. 0.1 * sqrt((p.rhog - p.rhoa) / p.rhoa * p.g * d50)
        Tr = @. (ustar_grain^2 - ustar_th_B^2) / ustar_th_B^2
        gamma_r = @. 1 + 1 / Tr
        return @. (d90 + 15 * gamma_r * d50 * Tr^1) / 30
    end
end

"""
    velocity_stress!(s, m)

Shear stress from shear velocity: `tau = rhoa * ustar^2`, with components along the shear
velocity.
"""
function velocity_stress!(s::AeolisState, m::AeolisModel)
    rhoa = m.p.rhoa
    @inbounds for k in eachindex(s.ustar)
        u = s.ustar[k]
        tau = rhoa * u^2
        if u > 0
            s.taus[k] = tau * s.ustars[k] / u
            s.taun[k] = tau * s.ustarn[k] / u
        end
        s.tau[k] = hypot(s.taus[k], s.taun[k])
        if u == 0
            s.taus[k] = 0
            s.taun[k] = 0
            s.tau[k] = 0
        end
    end
    return s
end

"""
    stress_velocity!(s, m)

Shear velocity from shear stress: `ustar = sqrt(tau / rhoa)`, with components along the
stress.
"""
function stress_velocity!(s::AeolisState, m::AeolisModel)
    rhoa = m.p.rhoa
    @inbounds for k in eachindex(s.tau)
        tau = s.tau[k]
        s.ustar[k] = sqrt(tau / rhoa)
        if tau > 0
            s.ustars[k] = s.ustar[k] * s.taus[k] / tau
            s.ustarn[k] = s.ustar[k] * s.taun[k] / tau
        end
        if tau == 0
            s.ustar[k] = 0
            s.ustars[k] = 0
            s.ustarn[k] = 0
        end
    end
    return s
end

"""
    compute_shear!(s, m)

Perturb the shear stress by the topography (AeoLiS `wind.shear`): the spectral method of
Weng et al. (1991) with a separation bubble for `method_shear = :fft`, or the rotated-grid
methods `:quasi2d` (Kroy et al. 2002) and `:duna2d` along wind-aligned transects.
"""
function compute_shear!(s::AeolisState, m::AeolisModel)
    p = m.p
    (p.process_shear && m.ny > 0 && m.shear !== nothing) || return s
    if p.method_shear == :fft
        sh = m.shear
        taus, taun, hsep = wind_shear!(sh, s.x, s.y, s.zb, s.taus, s.taun, s.uw[1, 1],
                                       s.udir[1, 1], p.process_separation, p.c_b, p.mu_b,
                                       s.taus0[1, 1], s.taun0[1, 1],
                                       p.sep_filter_iterations,
                                       p.zsep_y_filter)
        s.taus .= taus
        s.taun .= taun
        s.tau .= hypot.(s.taus, s.taun)
        stress_velocity!(s, m)
        if p.process_separation
            s.hsep .= hsep
            s.zsep .= s.hsep .+ s.zb
        end
    else
        tau, zsep, hsep = rotated_shear!(m.shear, s, m)
        s.tau .= tau
        @. s.taus = -s.tau * sin((-p.alfa + s.udir) / 180.0 * pi)
        @. s.taun = -s.tau * cos((-p.alfa + s.udir) / 180.0 * pi)
        stress_velocity!(s, m)
    end
    return s
end
