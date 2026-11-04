# Shear velocity threshold (AeoLiS `threshold.py`): grain size (Bagnold), moisture (several
# formulations), salt, sheltering by roughness elements and a non-erodible layer.

const MOISTURE_THRESHOLDS = (:belly_johnson, :hotta, :chepil, :saleh_fryear,
                             :saleh_fryear_mod, :shao, :dong_2002, :gregory_darwish,
                             :cornelis, :dong_2007)

"""
    compute_threshold!(s, m)

The shear velocity threshold of every cell and fraction from the enabled `th_*` processes,
the threshold mask, and the non-erodible layer. With `threshold_file` the threshold is fixed
(only the non-erodible layer applies).
"""
function compute_threshold!(s::AeolisState, m::AeolisModel)
    p = m.p
    if p.process_threshold && p.threshold_file === nothing
        p.th_grainsize && threshold_grainsize!(s, m)
        # th_bedslope has no effect in AeoLiS
        if p.th_moisture
            threshold_moisture!(s, m)
        elseif p.process_tide
            for j in axes(s.zb, 2), i in axes(s.zb, 1)
                if s.zb[i, j] - s.zs[i, j] < -p.eps
                    s.uth[i, j, :] .= Inf
                end
            end
        end
        p.th_salt && threshold_salt!(s, m)
        p.th_sheltering && threshold_sheltering!(s, m)
        for q in 1:m.nf
            @views apply_mask!(s.uth[:, :, q], m.masks[:threshold_mask])
        end
        s.uthf .= s.uth
    end
    p.th_nelayer && threshold_nelayer!(s, m)
    return s
end

"""
    threshold_grainsize!(s, m)

Bagnold's threshold `Aa * sqrt((rhog - rhoa) / rhoa * g * d)` for every fraction.
"""
function threshold_grainsize!(s::AeolisState, m::AeolisModel)
    p = m.p
    for q in 1:m.nf
        s.uth[:, :, q] .= p.Aa * sqrt((p.rhog - p.rhoa) / p.rhoa * p.g * m.grain_size[q])
    end
    s.uth0 .= s.uth
    return s
end

"""
    threshold_moisture!(s, m)

Raise the threshold with the gravimetric moisture content by `method_moist_threshold`, and
make it infinite above `max_moist` percent volumetric moisture.
"""
function threshold_moisture!(s::AeolisState, m::AeolisModel)
    p = m.p
    meth = _check_option(:method_moist_threshold, p.method_moist_threshold,
                         MOISTURE_THRESHOLDS)
    conv = p.rhow / (p.rhog * (1.0 - p.porosity))
    lim = p.resw_moist * conv + 0.005
    d = sum(p.grain_size .* vec(m.grain_dist[1, :]))
    K = if meth == :dong_2002
        dmm = d * 1000
        dmm < 0.135 ? 2.51 :
        dmm < 0.150 ? 2.05 :
        dmm < 0.200 ? 2.75 :
        dmm < 0.250 ? 1.59 : dmm < 0.400 ? 1.87 : 2.15
    else
        0.0
    end
    th_mg = p.max_moist / 100 * conv
    for q in 1:m.nf, j in axes(s.zb, 2), i in axes(s.zb, 1)
        mg = s.moist[i, j] * conv
        u = s.uth[i, j, q]
        if mg > lim
            u = if meth == :belly_johnson
                u * max(1.0, 1.8 + 0.6 * log10(mg * 100.0))
            elseif meth == :hotta
                u + 7.5 * mg
            elseif meth == :chepil
                sqrt(u^2 + 0.6 * mg^2 / (p.rhoa * p.w1_5^2))
            elseif meth == :saleh_fryear
                0.305 + 0.022 * mg / p.w1_5 + 0.506 * (mg / p.w1_5)^2
            elseif meth == :saleh_fryear_mod
                u + 0.022 * mg / p.w1_5 + 0.506 * (mg / p.w1_5)^2
            elseif meth == :shao
                u * exp(22.7 * mg)
            elseif meth == :dong_2002
                u * sqrt(1 + K * 100 * mg)
            elseif meth == :gregory_darwish
                wdiff = max(mg - 0.5 * p.w1_5, 0.0)
                u * sqrt(1 + mg + 6 / pi * 6.12e-05 / (p.rhog * p.g * d^2) +
                         738.8 / (p.rhow * p.g * d) * exp(-0.1 * mg / p.w1_5) * wdiff)
            elseif meth == :cornelis
                sqrt(0.013 *
                     (1 + mg +
                      1.7e-4 / ((p.rhog - p.rhoa) * p.g * d^2) *
                      (1 + 3e14 * 0.075^2 * d / (10^9 * exp(-6.5 * mg / p.w1_5)))) *
                     (p.rhog - p.rhoa) / p.rhoa * p.g * d)
            else # :dong_2007
                0.16 * sqrt((p.rhog - p.rhoa) / p.rhoa * p.g * d) * (1 + 478.2 * mg^1.52)
            end
        end
        mg > th_mg && (u = Inf)
        s.uth[i, j, q] = u
    end
    return s
end

"""
    threshold_salt!(s, m)

Raise the threshold with the salt content of the top layer (Nickling and Ecclestone 1981).
"""
function threshold_salt!(s::AeolisState, m::AeolisModel)
    p = m.p
    for q in 1:m.nf, j in axes(s.zb, 2), i in axes(s.zb, 1)
        cs = p.csalt * (1.0 - s.salt[i, j, 1])
        s.uth[i, j, q] *= 1.03 * exp(0.1027 * 1e3 * cs)
    end
    return s
end

"""
    threshold_sheltering!(s, m)

Raise the threshold by sheltering from non-erodible roughness elements (Raupach et al. 1993),
from the share of top-layer mass in fractions that the wind cannot move.
"""
function threshold_sheltering!(s::AeolisState, m::AeolisModel)
    p = m.p
    for j in axes(s.zb, 2), i in axes(s.zb, 1)
        gd = 0.0
        for q in 1:m.nf
            s.ustar[i, j] <= s.uth[i, j, q] && (gd += s.mass[i, j, 1, q])
        end
        gd /= sum(@view s.mass[i, j, 1, :])
        Rti = sqrt((1.0 - p.m * gd) * (1.0 + p.m * p.beta / p.sigma * gd))
        s.Rti[i, j] = Rti
        for q in 1:m.nf
            s.uth[i, j, q] *= Rti
        end
    end
    return s
end

"""
    threshold_nelayer!(s, m)

Make the threshold infinite where the bed is within 0.01 m of the non-erodible layer.
"""
function threshold_nelayer!(s::AeolisState, m::AeolisModel)
    s.zne .= m.ne
    for j in axes(s.zb, 2), i in axes(s.zb, 1)
        if s.zb[i, j] <= s.zne[i, j] + 0.01
            s.uth[i, j, :] .= Inf
        end
    end
    return s
end
