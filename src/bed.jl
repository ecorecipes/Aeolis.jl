# Bed composition and bed level (AeoLiS `bed.py`): initial layers and fractions, mixing of the
# top layer by waves, the bed update from pickup with the layer bookkeeping, resetting the
# intertidal bed, and the averaged bed level change used by vegetation.

const SECONDS_PER_YEAR = 3600.0 * 24.0 * 365.25

"""
    initialize_bed!(s, m)

Bed level, non-erodible layer, layer thickness and the sediment mass per layer and fraction
from `grain_dist` (one distribution, or one row per layer), `bedcomp_file` or `supply_file`.
"""
function initialize_bed!(s::AeolisState, m::AeolisModel)
    p = m.p
    nl, nf = m.nl, m.nf
    s.zb .= m.bed
    s.zb0 .= m.bed
    s.zne .= m.ne
    s.zdry .= 0.05
    s.thlyr .= p.layer_thickness
    if p.bedcomp_file !== nothing && p.supply_file !== nothing
        throw(ConfigurationError("bedcomp_file and supply_file cannot both be set"))
    end
    if p.supply_file !== nothing
        s.mass .= 0
    elseif m.bedcomp === nothing
        gd = m.grain_dist
        if size(gd, 1) > 1 && size(gd, 1) != nl
            throw(ConfigurationError("grain_dist has $(size(gd, 1)) rows; one per layer ($nl) is needed"))
        end
        for i in 1:nl
            gs = gd[size(gd, 1) == 1 ? 1 : i, :]
            gs = gs ./ sum(gs)
            for j in 1:nf
                @views s.mass[:, :, i, j] .= p.rhog * (1.0 - p.porosity) .*
                                             s.thlyr[:, :, i] .* gs[j]
            end
        end
    else
        s.mass .= m.bedcomp
    end
    if p.threshold_file !== nothing
        th = _as_grid(eltype(s.zb), p.threshold_file, m.nx + 1, m.ny + 1, :threshold_file)
        for j in 1:nf
            s.uth[:, :, j] .= th
        end
    end
    return s
end

"""
    mixtoplayer!(s, m)

Mix the layers within the depth of disturbance (`facDOD` times the wave height for mixing)
when `process_mixtoplayer` is on.
"""
function mixtoplayer!(s::AeolisState, m::AeolisModel)
    p = m.p
    p.process_mixtoplayer || return s
    nl, nf = m.nl, m.nf
    DOD = p.facDOD .* s.Hsmix
    th = dropdims(sum(s.thlyr; dims=3); dims=3)
    f = ones(size(DOD))
    for k in eachindex(DOD)
        DOD[k] > 0 && (f[k] = min(1.0, th[k] / DOD[k]))
    end
    cth = cumsum(s.thlyr; dims=3)
    ix = falses(size(s.mass))
    for l in 1:nl, j in axes(DOD, 2), i in axes(DOD, 1), q in 1:nf
        ix[i, j, l, q] = (cth[i, j, l] <= DOD[i, j]) && (DOD[i, j] > 0)
    end
    if any(ix)
        ix[:, :, 1, :] .= true
        for q in 1:nf, j in axes(DOD, 2), i in axes(DOD, 1)
            acc, cnt = 0.0, 0
            for l in 1:nl
                if ix[i, j, l, q]
                    acc += s.mass[i, j, l, q]
                    cnt += 1
                end
            end
            mean1 = acc / cnt
            for l in 1:nl
                if ix[i, j, l, q]
                    s.mass[i, j, l, q] = mean1 * f[i, j] +
                                         s.mass[i, j, l, q] * (1 - f[i, j])
                end
            end
        end
    end
    return s
end

"""
    wet_bed_reset!(s, m)

Relax the bed below the total water level towards its initial level with time scale
`Tbedreset` when `process_wet_bed_reset` is on.
"""
function wet_bed_reset!(s::AeolisState, m::AeolisModel)
    m.p.process_wet_bed_reset || return s
    r = m.dt_opt / m.p.Tbedreset
    @inbounds for k in eachindex(s.zb)
        if s.TWL[k] > s.zb[k]
            s.zb[k] += (s.zb0[k] - s.zb[k]) * r
        end
    end
    return s
end

"""
    update_bed!(s, m)

Remove the pickup from the top layer, shift mass between layers so that layer thicknesses are
kept, take new sediment from below the bottom layer with the initial composition, and change
the bed level by the net deposition (AeoLiS `bed.update`). Like AeoLiS this may modify
`s.pickup` where it would remove more sediment than available.
"""
function update_bed!(s::AeolisState, m::AeolisModel)
    p = m.p
    nx, ny, nl, nf = m.nx, m.ny, m.nl, m.nf
    if m.supply !== nothing
        @views s.mass[:, :, 1, 1] .-= s.pickup[:, :, 1]
        @views s.mass[:, :, 1, 1] .+= m.supply .* m.dt_opt
        if p.process_tide
            for k in CartesianIndices(s.zb)
                s.zb[k] < s.zs[k] && (s.mass[k, 1, 1] = 0)
            end
        end
        return s
    end
    ncell = (nx + 1) * (ny + 1)
    pickup = reshape(s.pickup, ncell, nf)            # a view: changes reach s.pickup
    dm = repeat(-sum(pickup; dims=2), 1, nf)
    ix_ero = dm[:, 1] .< 0
    ix_dep = dm[:, 1] .> 0
    mm = copy(reshape(s.mass, ncell, nl, nf))
    _prevent_negative_mass!(mm, dm, pickup)
    d = normalize_along(mm, 3)
    @views mm[:, 1, :] .-= pickup
    _arrange_layers!(mm, dm, d, nl, ix_ero, ix_dep)
    gd = m.grain_dist
    if size(gd, 1) > 1
        g = gd[end, :] ./ sum(gd[end, :])
        for c in findall(ix_ero), q in 1:nf
            mm[c, nl, q] -= dm[c, q] * g[q]
        end
    elseif m.bedcomp !== nothing
        bc = reshape(m.bedcomp, ncell, nl, nf)
        for c in findall(ix_ero)
            b = bc[c, nl, :]
            sb = sum(b)
            for q in 1:nf
                mm[c, nl, q] -= dm[c, q] * (sb != 0 ? b[q] / sb : 0.0)
            end
        end
    else
        g = gd[1, :] ./ sum(gd[1, :])
        for c in findall(ix_ero), q in 1:nf
            mm[c, nl, q] -= dm[c, q] * g[q]
        end
    end
    prevent_tiny_negatives!(mm, p.max_error)
    if minimum(mm) < 0
        @warn "negative mass" nrcells = count(c -> any(<(0), mm[c, :, :]), 1:ncell) minvalue = minimum(mm) time = m.time
    end
    s.mass .= reshape(mm, size(s.mass))
    @views s.masstop .= s.mass[:, :, 1, :]
    if p.process_bedupdate
        dz = reshape(dm[:, 1], nx + 1, ny + 1) ./ (p.rhog * (1.0 - p.porosity))
        s.dzb .= dz
        s.zb .+= dz
        p.process_tide && (s.zs .+= dz)
    end
    return s
end

# AeoLiS `prevent_negative_mass`: redistribute erosion and deposition within a cell, and
# shift layers when deposition exceeds the top layer. `pickup` and `dm` change in place.
function _prevent_negative_mass!(mm, dm, pickup)
    ncell, nl, nf = size(mm)
    ix_dep = dm[:, 1] .> 0
    ero = max.(0.0, pickup)
    dep = -min.(0.0, pickup)
    erog = sum(ero; dims=2)
    ix = ix_dep .& (erog[:, 1] .> 0)
    if any(ix)
        d = normalize_along(dep, 2)
        for c in findall(ix)
            ddep = erog[c, 1] .* d[c, :]
            pickup[c, :] .= -dep[c, :] .+ ddep
            dm[c, :] .= -sum(pickup[c, :])
            mm[c, 1, :] .-= ero[c, :] .- ddep
        end
    end
    mx = sum(mm[:, 1, :]; dims=2)
    n = [_floordiv(dm[c, 1], mx[c, 1]) for c in 1:ncell]
    if any(>(0), n)
        d = normalize_along(pickup, 2)
        for i in 1:nl
            ixs = findall(>(i - 1), n)
            isempty(ixs) && break
            for c in ixs
                mm[c, (i + 1):nl, :] .= mm[c, i:(nl - 1), :]
                mm[c, i, :] .= mx[c, 1] .* d[c, :]
                pickup[c, :] .-= mm[c, i, :]
            end
        end
        ix2 = dm[:, 1] .> mx[:, 1]
        for c in findall(ix2)
            pickup[c, :] .= 0
            dm[c, :] .= -sum(pickup[c, :])
        end
    end
    return mm
end

# numpy floor division of floats, including division by zero
function _floordiv(a::AbstractFloat, b::AbstractFloat)
    if b == 0
        return a == 0 || isnan(a) ? NaN : copysign(Inf, a)
    end
    return floor(a / b)
end

function _arrange_layers!(mm, dm, d, nl, ix_ero, ix_dep)
    ero = findall(ix_ero)
    dep = findall(ix_dep)
    nf = size(mm, 3)
    for i in 2:nl
        for c in ero, q in 1:nf
            mm[c, i - 1, q] -= dm[c, q] * d[c, i, q]
            mm[c, i, q] += dm[c, q] * d[c, i, q]
        end
        for c in dep, q in 1:nf
            mm[c, i - 1, q] -= dm[c, q] * d[c, i - 1, q]
            mm[c, i, q] += dm[c, q] * d[c, i - 1, q]
        end
    end
    for c in dep, q in 1:nf
        mm[c, nl, q] -= dm[c, q] * d[c, nl, q]
    end
    return mm
end

"""
    average_change!(s, m)

Bed level change over the step, per year, and averaged over `avg_time` (AeoLiS
`average_change`); `dzbveg` is zero until `avg_time` has passed.
"""
function average_change!(s::AeolisState, m::AeolisModel)
    p = m.p
    l = m.prev
    s.dzb .= s.zb .- l.zb
    s.dzbyear .= s.dzb .* SECONDS_PER_YEAR ./ (m.dt_opt * p.accfac)
    n = (m.dt_opt * p.accfac) / p.avg_time
    s.dzbavg .= n .* s.dzbyear .+ (1 - n) .* l.dzbavg
    s.dzbveg .= s.dzbavg
    m.time < p.avg_time && (s.dzbveg .= 0)
    return s
end
