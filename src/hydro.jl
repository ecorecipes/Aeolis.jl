# Hydrodynamics and moisture (AeoLiS `hydro.py`): water levels from tides, wave heights and
# runup (Stockdon et al. 2006), the groundwater level (Boussinesq equation), surface moisture
# by infiltration or by soil water retention with hysteresis, and salt.

"""
    pyrange(n, start, stop)

The 1-based indices of the Python slice `start:stop` of a sequence of length `n`, with
negative bounds counting from the end, as numpy slices them.
"""
function pyrange(n::Int, start::Int, stop::Int)
    a = start < 0 ? max(start + n, 0) : min(start, n)
    b = stop < 0 ? max(stop + n, 0) : min(stop, n)
    return (a + 1):b
end

"""
    stockdon_runup(Ho, Tp, beta)

Setup, swash and runup `(eta, sigma_s, R)` of waves of height `Ho` and period `Tp` on a beach
of slope `beta` (Stockdon et al. 2006, as in AeoLiS).
"""
function stockdon_runup(Ho::Real, Tp::Real, beta::Real)
    if Ho > 0 && Tp > 0 && beta > 0
        Lo = 9.81 * Tp * Tp / (2 * pi)
        iribarren = beta / (Ho / Lo)^0.5
        if iribarren < 0.3
            R = 0.043 * sqrt(Ho * Lo)
            sigma_s = 0.046 * sqrt(Ho * Lo) / 2
            eta = R - sigma_s
        else
            nsigma = 2
            eta = 0.35 * beta * sqrt(Ho * Lo)
            sigma_s = sqrt(Ho * Lo * (0.563 * (beta * beta) + 0.0004)) * nsigma / 2 / 2
            R = 1.1 * (eta + sigma_s)
        end
        return eta, sigma_s, R
    end
    return 0.0, 0.0, 0.0
end

# Grid version (AeoLiS applies it to externally supplied wave fields)
function stockdon_runup(Ho::AbstractArray, Tp::AbstractArray, beta::Real)
    R = zeros(size(Ho))
    sigma_s = zeros(size(Ho))
    eta = zeros(size(Ho))
    for k in eachindex(Ho)
        Ho[k] > 0 || continue
        Lo = 9.81 * Tp[k] * Tp[k] / (2 * pi)
        iri = beta / (Ho[k] / Lo)^0.5
        if iri < 0.3
            R[k] = 0.043 * sqrt(Ho[k] * Lo)
            sigma_s[k] = 0.046 * sqrt(Ho[k] * Lo) / 2
            eta[k] = R[k] - sigma_s[k]
        elseif iri > 0.3
            eta[k] = 0.35 * beta * sqrt(Ho[k] * Lo)
            sigma_s[k] = sqrt(Ho[k] * Lo * (0.563 * beta^2 + 0.0004)) * 2 / 2 / 2
            R[k] = 1.1 * (eta[k] + sigma_s[k])
        end
    end
    return eta, sigma_s, R
end

"""
    interpolate_hydro!(s, m, t)

Still and total water levels, wave height and period, the wave height for mixing, and the
meteorological conditions at time `t`; the water level `zs` is at least the bed.
"""
function interpolate_hydro!(s::AeolisState, m::AeolisModel, t::Real)
    p = m.p
    if p.process_tide
        if :zs ∉ p.external_vars
            s.SWL .= m.tide === nothing ? 0.0 : m.tide(t)
            apply_mask!(s.SWL, m.masks[:tide_mask])
        else
            s.SWL .= s.zs
            apply_mask!(s.SWL, m.masks[:tide_mask])
            _nearest_fill!(s.SWL, (s.zs .- s.zb) .> 2 * p.eps)
            @warn "SWL is set equal to the external zs (DSWL = zs)" maxlog = 1
        end
    else
        s.SWL .= s.zb .* 0
    end
    if :Hs ∉ p.external_vars && :Tp ∉ p.external_vars
        if p.process_wave && m.wave_height !== nothing
            h = max.(0.0, s.SWL .- s.zb)
            s.Hs .= m.wave_height(t)
            s.Tp .= m.wave_period(t)
            if p.process_runup
                rm = m.masks[:runup_mask]
                for iy in axes(s.zb, 2)
                    eta, sigma_s, R = stockdon_runup(s.Hs[1, iy], s.Tp[1, iy],
                                                     p.beach_slope)
                    s.R[:, iy] .= R
                    s.eta[:, iy] .= eta
                    s.sigma_s[:, iy] .= sigma_s
                    if rm isa AbstractArray
                        @views apply_mask!(s.eta[:, iy], rm[:, iy])
                        @views apply_mask!(s.R[:, iy], rm[:, iy])
                    end
                    @views s.TWL[:, iy] .= s.SWL[:, iy] .+ s.R[:, iy]
                    @views s.DSWL[:, iy] .= s.SWL[:, iy] .+ s.eta[:, iy]
                end
            end
            s.Hs .= min.(h .* p.gamma, s.Hs)
            apply_mask!(s.Hs, m.masks[:wave_mask])
            apply_mask!(s.Tp, m.masks[:wave_mask])
        else
            s.Hs .= s.zb .* 0
            s.Tp .= s.zb .* 0
        end
    else
        apply_mask!(s.Hs, m.masks[:wave_mask])
        apply_mask!(s.Tp, m.masks[:wave_mask])
    end
    if p.process_runup && :Hs in p.external_vars
        eta, sigma_s, R = stockdon_runup(s.Hs, s.Tp, p.beach_slope)
        s.R .= R
        rm = m.masks[:runup_mask]
        if rm isa AbstractArray
            apply_mask!(s.eta, rm)
            apply_mask!(s.R, rm)
        end
        s.TWL .= s.SWL .+ s.R
        s.DSWL .= s.SWL
    end
    if p.process_wave && m.wave_height !== nothing
        h_mix = max.(0.0, s.TWL .- s.zb)
        s.Hsmix .= m.wave_height(t)
        s.Hsmix .= min.(h_mix .* p.gamma, s.Hsmix)
        apply_mask!(s.Hsmix, m.masks[:wave_mask])
    end
    if p.process_moist && p.method_moist_process == :surf_moisture && m.meteo !== nothing
        vals = [f(t) for f in m.meteo]
        for (k, name) in enumerate((:T, :Q, :RH, :P, :U))
            m.meteo_values[name] = vals[k]
        end
    end
    s.zs .= s.SWL
    for k in eachindex(s.zs)
        s.zb[k] > s.zs[k] && (s.zs[k] = s.zb[k])
    end
    return s
end

# Replace the values outside `valid` by the nearest valid value (scipy NearestNDInterpolator
# on grid indices).
function _nearest_fill!(a::AbstractMatrix, valid::AbstractMatrix)
    pts = findall(valid)
    isempty(pts) && return a
    for k in CartesianIndices(a)
        valid[k] && continue
        best = argmin(q -> (q[1] - k[1])^2 + (q[2] - k[2])^2, pts)
        a[k] = a[best]
    end
    return a
end

"""
    update_hydro!(s, m, dt, t)

Groundwater (fourth-order Runge-Kutta on the Boussinesq equation in sub-steps of
`dt / tfac_gw`, with runup overheight), moisture (`:infiltration` or `:surf_moisture` with
soil water retention, optionally with scanning curves), and salt.
"""
function update_hydro!(s::AeolisState, m::AeolisModel, dt::Real, t::Real)
    p = m.p
    if p.process_groundwater
        if t == p.tstart
            s.gw .= p.in_gw
            s.gw_prev .= s.gw
        end
        s.wetting .= s.gw .> s.gw_prev
        s.gw_prev .= s.gw      # AeoLiS aliases gw_prev to gw; see the first sub-step below
        dt_gw = trunc(Int, dt / p.tfac_gw)
        bound = _check_option(:boundary_gw, p.boundary_gw, (:no_flow, :static)) ==
                :no_flow ?
                0 : 1
        for i in 0:(trunc(Int, dt / dt_gw) - 1)
            t_gw = t + i * dt_gw
            interpolate_hydro!(s, m, t_gw)
            shl_ix = [something(findfirst(view(s.zb, :, r) .> view(s.DSWL, :, r)), 1) - 2
                      for r in axes(s.zb, 2)]
            runup_ix = [something(findfirst(view(s.zb, :, r) .> view(s.TWL, :, r)), 1) - 2
                        for r in axes(s.zb, 2)]
            f1 = _boussinesq!(s.gw, s.DSWL, s.ds, p.GW_stat, p.K_gw, p.ne_gw, p.D_gw,
                              shl_ix,
                              bound, s.zb, p.process_seepage_face)
            # the first call changes gw in place, and in AeoLiS gw_prev is the same array
            i == 0 && (s.gw_prev .= s.gw)
            f2 = _boussinesq!(s.gw .+ dt_gw / 2 .* f1, s.DSWL, s.ds, p.GW_stat, p.K_gw,
                              p.ne_gw, p.D_gw, shl_ix, bound, s.zb, p.process_seepage_face)
            f3 = _boussinesq!(s.gw .+ dt_gw / 2 .* f2, s.DSWL, s.ds, p.GW_stat, p.K_gw,
                              p.ne_gw, p.D_gw, shl_ix, bound, s.zb, p.process_seepage_face)
            f4 = _boussinesq!(s.gw .+ dt_gw .* f3, s.DSWL, s.ds, p.GW_stat, p.K_gw, p.ne_gw,
                              p.D_gw, shl_ix, bound, s.zb, p.process_seepage_face)
            s.gw .= s.gw .+ dt_gw / 6 .* (f1 .+ 2 .* f2 .+ 2 .* f3 .+ f4)
            if p.process_wave
                fx = zeros(size(s.gw))
                _runup_overheight!(fx, shl_ix, runup_ix, s.x)
                s.gw .+= p.Cl_gw .* fx
            end
            apply_mask!(s.gw, m.masks[:gw_mask])
            s.gw .= min.(s.gw, s.zb)
            for k in eachindex(s.gw)
                s.zb[k] < s.DSWL[k] && (s.gw[k] = s.DSWL[k])
            end
        end
    end
    wet = nothing
    if p.process_moist
        meth = _check_option(:method_moist_process, p.method_moist_process,
                             (:infiltration, :surf_moisture))
        F1 = -log(0.5) / p.Tdry
        if meth == :infiltration
            wet = (s.TWL .- s.zb) .> p.eps
            for k in eachindex(s.moist)
                s.moist[k] = wet[k] ? p.porosity : s.moist[k] * exp(-F1 * dt)
                s.moist[k] = max(0.0, min(p.porosity, s.moist[k]))
            end
        else
            _surf_moisture!(s, m, dt, t, F1)
            wet = s.TWL .> s.zb
        end
    end
    if p.process_salt
        wet === nothing &&
            throw(ConfigurationError("process_salt needs process_moist, which defines the wet cells"))
        met = m.meteo_values
        F2 = -log(0.5) / p.Tsalt
        pcp = met[:RH] / 3600.0 / 1000.0
        for j in axes(s.zb, 2), i in axes(s.zb, 1)
            s.salt[i, j, 1] = wet[i, j] ? 1.0 : s.salt[i, j, 1] * exp(-F2 * dt)
            s.salt[i, j, 1] = min(1.0, s.salt[i, j, 1] + pcp * dt / p.layer_thickness)
        end
    end
    return s
end

function _surf_moisture!(s, m, dt, t, F1)
    p = m.p
    s.moist .+= min.(0.0, (s.moist .- p.fc) .* (exp(-F1 * dt) - 1))
    wet = s.TWL .> s.zb
    s.moist[wet] .= p.satd_moist
    met = m.meteo_values
    evo = _evaporation(s, p, met) ./ 24.0 ./ 3600.0 ./ 1000.0
    pcp = met[:RH] / 3600.0 / 1000.0
    for k in eachindex(s.moist)
        wet[k] && continue
        s.moist[k] = max(s.moist[k] + (pcp - evo[k]) * dt / p.thick_moist, p.resw_moist)
        s.moist[k] = min(s.moist[k], p.satd_moist)
    end
    h = max.(0.0, (s.zb .- s.gw) .* 100)
    vg(res, sat, alfa, n, mm, hh) = res + (sat - res) / (1 + abs(alfa * hh)^n)^mm
    if p.process_scanning
        if t == 0
            s.moist_swr .= vg.(p.resw_moist, p.satw_moist, p.alfaw_moist, p.nw_moist,
                               p.mw_moist, h)
            s.h_delta .= max.(0.0, (s.zb .- s.gw) .* 100)
        else
            _hdelta!(s, p)
            s.w_h .= vg.(p.resw_moist, p.satw_moist, p.alfaw_moist, p.nw_moist, p.mw_moist,
                         h)
            s.w_hdelta .= vg.(p.resw_moist, p.satw_moist, p.alfaw_moist, p.nw_moist,
                              p.mw_moist, s.h_delta)
            s.d_h .= vg.(p.resd_moist, p.satd_moist, p.alfad_moist, p.nd_moist, p.md_moist,
                         h)
            s.d_hdelta .= vg.(p.resd_moist, p.satd_moist, p.alfad_moist, p.nd_moist,
                              p.md_moist, s.h_delta)
            @. s.scan_w_moist = max(min(s.w_h +
                                        (p.satw_moist - s.w_h) /
                                        max(p.satw_moist - s.w_hdelta, 0.0001) *
                                        (s.d_hdelta - s.w_hdelta), s.d_h), s.w_h)
            @. s.scan_d_moist = max(min(s.w_h +
                                        (s.w_hdelta - s.w_h) /
                                        max(p.satd_moist - s.w_h, 0.0001) *
                                        (s.d_h - s.w_h), s.d_h), s.w_h)
            _swr_curve!(s)
        end
    else
        for k in eachindex(h)
            s.moist_swr[k] = s.wetting[k] > 0 ?
                             vg(p.resw_moist, p.satw_moist, p.alfaw_moist, p.nw_moist,
                                p.mw_moist, h[k]) :
                             vg(p.resd_moist, p.satd_moist, p.alfad_moist, p.nd_moist,
                                p.md_moist, h[k])
        end
    end
    s.moist .= min.(max.(s.moist, s.moist_swr), p.satd_moist)
    return s
end

function _evaporation(s, p, met)
    l = 2.26
    sp(T) = 0.6108 * exp(17.27 * T / (T + 237.3))
    mslope = 4098.0 * sp(met[:T]) / (met[:T] + 237.3)^2
    delta = sp(met[:T]) * (1.0 - met[:U] / 100)
    gamma = (p.cpair * met[:P]) / (0.622 * l)
    u2 = 0.174 / log10(p.z / 2.0) .* s.uw
    return @. (mslope * met[:Q] + 6.43 * gamma * delta * (1.0 + 0.86 * u2)) /
              (l * (mslope + gamma))
end

# AeoLiS `Boussinesq`: the rate of change of the groundwater level landward of the shoreline.
# Changes GW in place at the shoreline and the landward boundary, as AeoLiS does. Indices are
# 0-based as in AeoLiS; `pyrange` turns slices into ranges.
function _boussinesq!(GW, DSWL, ds, GW_stat, K, ne, D, shl_ix, bound, zb, seepage)
    n, nrow = size(GW)
    for r in 1:nrow
        a = _py(shl_ix[r], n) + 1
        GW[a, r] = DSWL[a, r]
        rp = r == 1 ? nrow : r - 1
        b = _py(shl_ix[rp], n) + 1
        GW[b, r] = DSWL[b, r]
    end
    if bound == 0
        for r in 1:nrow
            v = GW[n - 3, r]
            GW[n, r] = v
            GW[n - 1, r] = v
            GW[n - 2, r] = v
        end
    else
        GW[n, :] .= GW_stat
        GW[n - 1, :] .= GW_stat
        GW[n - 2, :] .= GW_stat
    end
    if seepage
        for r in 1:nrow
            ixs = something(findfirst(k -> !(GW[k, r] + 0.05 >= zb[k, r]), 1:n), 1) - 1
            if shl_ix[r] < ixs - 1
                rr = pyrange(n, shl_ix[r], ixs - 1)
                GW[rr, r] .= zb[rr, r]
            end
        end
    end
    dGW = zeros(size(GW))
    a = zeros(size(GW))
    b = zeros(size(GW))
    c = zeros(size(GW))
    for r in 1:nrow
        sh = shl_ix[r]
        sh < n - 3 || continue
        c0 = pyrange(n, sh, -2)
        @views a[c0, r] .= (GW[pyrange(n, sh + 1, -1), r] .- 2 .* GW[c0, r] .+
                            GW[pyrange(n, sh - 1, -3), r]) ./ ds[c0, r] .^ 2
        @views b[c0, r] .= (GW[c0, r] .* (GW[pyrange(n, sh + 1, -1), r] .+
                                          GW[pyrange(n, sh - 1, -3), r])) ./ ds[c0, r]
        c1 = pyrange(n, sh + 1, -3)
        @views c[c1, r] .= (b[pyrange(n, sh + 2, -2), r] .- b[pyrange(n, sh, -4), r]) ./
                           ds[c1, r]
        @views dGW[c1, r] .= K / ne .* (D .* a[c1, r] .+ c[c1, r])
    end
    return dGW
end

function _runup_overheight!(fx, shl_ix, runup_ix, x)
    n, nrow = size(fx)
    for r in 1:nrow
        sh, ru = shl_ix[r], runup_ix[r]
        fxi = trunc(Int, sh + 2 / 3 * (ru - sh))
        xs = x[_py(sh, n) + 1, r]
        xr = x[_py(ru, n) + 1, r]
        for k in pyrange(n, sh, fxi)
            fx[k, r] = (x[k, r] - xs) / (2 / 3 * (xr - xs))
        end
        for k in pyrange(n, fxi + 1, ru)
            fx[k, r] = 3 - (x[k, r] - xs) / (1 / 3 * (xr - xs))
        end
        fx[_py(fxi, n) + 1, r] = 1
    end
    return fx
end

function _hdelta!(s, p)
    for k in eachindex(s.h_delta)
        wetting = s.wetting[k] > 0
        scan_w = s.scan_w[k] > 0
        scan_d = s.scan_d[k] > 0
        gw, gwp = s.gw[k], s.gw_prev[k]
        if (!scan_w && wetting && gw < gwp) || (!scan_d && !wetting && gw > gwp)
            s.h_delta[k] = max(0.0, (s.zb[k] - gw) * 100)
        end
        if scan_w && wetting && gw < gwp
            w_int = min((s.scan_w_moist[k] - s.w_h[k]) * (p.satd_moist - s.w_h[k]) /
                        (s.d_h[k] - s.w_h[k]) + s.w_h[k], p.satw_moist)
            s.h_delta[k] = max(1 / p.alfaw_moist *
                               (((p.satw_moist - p.resw_moist) /
                                 max(w_int - p.resw_moist, 0.00001))^(1 / p.mw_moist) - 1)^(1 /
                                                                                            p.nw_moist),
                               0.0)
        end
        if scan_d && !wetting && gw > gwp
            hit = 0.0
            F = 1.0
            while F > 0.01
                hit += 0.01
                wd = p.resw_moist +
                     (p.satw_moist - p.resw_moist) /
                     (1 + abs(p.alfaw_moist * hit)^p.nw_moist)^p.mw_moist
                dd = p.resd_moist +
                     (p.satd_moist - p.resd_moist) /
                     (1 + abs(p.alfad_moist * hit)^p.nd_moist)^p.md_moist
                F = s.w_h[k] +
                    (p.satw_moist - s.w_h[k]) / max(p.satw_moist - wd, 0.0001) *
                    (dd - wd) - s.scan_d_moist[k]
            end
            s.h_delta[k] = hit
        end
    end
    return s
end

function _swr_curve!(s)
    for k in eachindex(s.moist_swr)
        gw, gwp = s.gw[k], s.gw_prev[k]
        wetting = s.wetting[k] > 0
        scan_w = s.scan_w[k] > 0
        scan_d = s.scan_d[k] > 0
        if gw >= gwp && wetting && !scan_w
            s.moist_swr[k] = s.w_h[k]
            s.scan_w[k] = 0
            s.scan_d[k] = 0
        elseif gw >= gwp && !wetting
            s.moist_swr[k] = s.scan_w_moist[k]
            s.scan_w[k] = s.scan_w_moist[k] > s.w_h[k]
            s.scan_d[k] = 0
        elseif gw >= gwp && wetting && scan_w
            s.moist_swr[k] = s.scan_w_moist[k]
            s.scan_w[k] = s.scan_w_moist[k] > s.w_h[k]
            s.scan_d[k] = 0
        elseif gw < gwp && !wetting && !scan_d
            s.moist_swr[k] = s.d_h[k]
            s.scan_d[k] = 0
            s.scan_w[k] = 0
        elseif gw < gwp && wetting
            s.moist_swr[k] = s.scan_d_moist[k]
            s.scan_d[k] = s.scan_d_moist[k] < s.d_h[k]
            s.scan_w[k] = 0
        elseif gw < gwp && !wetting && scan_d
            s.moist_swr[k] = s.scan_d_moist[k]
            s.scan_d[k] = s.scan_d_moist[k] < s.d_h[k]
            s.scan_w[k] = 0
        end
    end
    return s
end
