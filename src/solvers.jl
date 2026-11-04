# Transport solvers (AeoLiS `model.AeoLiS.solve*` and `utils.sweep`). Each advances the
# sediment concentration in the air `Ct` over one time step with an advection equation
# relaxing towards the saturated concentration, limited by the sediment available in the
# top layer of the bed, and returns the new concentration, fluxes, pickup and weights.
#
# The sparse matrices are assembled with `diags`, which places diagonals as
# `scipy.sparse.diags` does; the AeoLiS solvers take some diagonals from shifted copies of the
# coefficient arrays and that placement is part of the model.

# The weights at the first step: the initial grain size distribution (first row when there
# is one per layer), or the bed composition when it is given (AeoLiS uses `t == 0`).
function _initial_weights(s, m, w_init)
    if m.bedcomp !== nothing
        return copy(w_init)
    end
    g = m.grain_dist[1, :]
    return repeat(reshape(g, 1, 1, :), m.nx + 1, m.ny + 1, 1)
end

function _linsolve_cache(A, b, linsolve)
    prob = LinearSolve.LinearProblem(A, b)
    return linsolve === nothing ? LinearSolve.init(prob) : LinearSolve.init(prob, linsolve)
end

function _linsolve!(cache, b)
    cache.b = b
    sol = LinearSolve.solve!(cache)
    return copy(sol.u)
end

"""
    solve_trunk(s, m, dt, alpha, beta, linsolve)

The default AeoLiS solver: an implicit (`alpha = 1`), Crank-Nicolson (`alpha = 0.5`) or
explicit upwind discretisation solved as a sparse linear system per fraction, iterating the
fraction weights until the pickup is supported by the bed.
"""
function solve_trunk(s::AeolisState, m::AeolisModel, dt::Real, alpha::Real, beta::Real,
                     linsolve)
    p = m.p
    l = m.prev
    nxp, nyp = m.nx + 1, m.ny + 1
    nf = m.nf
    Ct = copy(s.Ct)
    pickup = copy(s.pickup)
    w_init, w_air, w_bed = compute_weights(s, m)
    w = m.time == 0 ? _initial_weights(s, m, w_init) : copy(w_init)
    Ti = dt / p.T
    for i in 1:nf
        us = s.us[:, :, i]
        un = s.un[:, :, i]
        us_plus = zeros(nxp, nyp)
        un_plus = zeros(nxp, nyp)
        us_min = zeros(nxp, nyp)
        un_min = zeros(nxp, nyp)
        us_plus[2:end, :] .= us[1:(end - 1), :]
        un_plus[:, 2:end] .= un[:, 1:(end - 1)]
        us_min[1:(end - 1), :] .= us[2:end, :]
        un_min[:, 1:(end - 1)] .= un[:, 2:end]
        us_plus[1, :] .= us[1, :]
        un_plus[:, 1] .= un[:, 1]
        us_min[end, :] .= us[end, :]
        un_min[:, end] .= un[:, end]
        Cs = dt .* s.dn .* s.dsdni .* us
        Cn = dt .* s.ds .* s.dsdni .* un
        Cs_plus = dt .* s.dn .* s.dsdni .* us_plus
        Cn_plus = dt .* s.ds .* s.dsdni .* un_plus
        Cs_min = dt .* s.dn .* s.dsdni .* us_min
        Cn_min = dt .* s.ds .* s.dsdni .* un_min
        b = abs(beta)
        if b >= 1
            ixs = Float64.(us .>= 0)
            ixn = Float64.(un .>= 0)
            sgs = 2 .* ixs .- 1
            sgn = 2 .* ixn .- 1
        else
            ixs = fill(b, nxp, nyp)
            ixn = fill(b, nxp, nyp)
            sgs = zeros(nxp, nyp)
            sgn = zeros(nxp, nyp)
        end
        A0 = 1 .+ (sgs .* Cs .+ sgn .* Cn .+ Ti) .* alpha
        Apx = Cn_min .* alpha .* (1 .- ixn)
        Ap1 = Cs_min .* alpha .* (1 .- ixs)
        Amx = -Cn_plus .* alpha .* ixn
        Am1 = -Cs_plus .* alpha .* ixs
        Ap2 = zeros(nxp, nyp)
        Am2 = zeros(nxp, nyp)
        A0[1, :] .= 1
        Apx[1, :] .= 0
        Amx[1, :] .= 0
        Am2[1, :] .= 0
        Am1[1, :] .= 0
        A0[end, :] .= 1
        Apx[end, :] .= 0
        Ap1[end, :] .= 0
        Ap2[end, :] .= 0
        Amx[end, :] .= 0
        bo = p.boundary_offshore
        if bo in (:flux, :noflux, :constant)
            Ap2[1, :] .= 0
            Ap1[1, :] .= 0
        elseif bo == :uniform
            Ap2[1, :] .= 0
            Ap1[1, :] .= -1
        elseif bo == :gradient
            Ap2[1, :] .= s.ds[2, :] ./ s.ds[3, :]
            Ap1[1, :] .= -1 .- s.ds[2, :] ./ s.ds[3, :]
        else
            throw(UnsupportedOptionError(:boundary_offshore, bo,
                                         "a cross-shore circular boundary is not implemented in the trunk solver"))
        end
        bn = p.boundary_onshore
        if bn == :flux || bo == :noflux   # AeoLiS tests the offshore setting here
            Am2[end, :] .= 0
            Am1[end, :] .= 0
        elseif bn == :constant
            Am2[end, :] .= 0
            Am1[end, :] .= 0
        elseif bn == :uniform
            Am2[end, :] .= 0
            Am1[end, :] .= -1
        elseif bn == :gradient
            Am2[end, :] .= s.ds[end - 1, :] ./ s.ds[end - 2, :]
            Am1[end, :] .= -1 .- s.ds[end - 1, :] ./ s.ds[end - 2, :]
        else
            throw(UnsupportedOptionError(:boundary_onshore, bn,
                                         "a cross-shore circular boundary is not implemented in the trunk solver"))
        end
        bl = p.boundary_lateral
        if bl == :constant
            for c in (1, nyp)
                A0[:, c] .= 1
                Apx[:, c] .= 0
                Ap1[:, c] .= 0
                Amx[:, c] .= 0
                Am1[:, c] .= 0
            end
        elseif bl != :circular
            throw(UnsupportedOptionError(:boundary_lateral, bl,
                                         "only circular and constant lateral boundaries are implemented in the trunk solver"))
        end
        n = nxp * nyp
        j = nxp
        ny = m.ny
        A = diags((vec(Apx)[1:j], vec(Amx)[(j + 1):end], vec(Am2)[3:end], vec(Am1)[2:end],
                   vec(A0), vec(Ap1)[1:(end - 1)], vec(Ap2)[1:(end - 2)],
                   vec(Apx)[(j + 1):end], vec(Amx)[1:j]),
                  (-j * ny, -j, -2, -1, 0, 1, 2, j, j * ny), n)
        renormalize_weights!(w, i)
        cache = nothing
        Ct_i = zeros(n)
        pickup_i = zeros(n)
        conv = false
        ixc = falses(n)
        mass_i = vec(s.mass[:, :, 1, i])
        Cu_i = vec(s.Cu[:, :, i])
        for it in 1:(p.max_iter)
            _count!(m, :matrixsolve)
            S = zeros(nxp, nyp)
            for q in 1:nf, jj in 1:nyp, ii in 1:nxp
                s.Cu[ii, jj, q] > 0 && (S[ii, jj] += s.Ct[ii, jj, q] / s.Cu[ii, jj, q])
            end
            s.S .= S
            I = 2:(nxp - 1)
            y_im = zeros(nxp, nyp)
            y_ex = zeros(nxp, nyp)
            @views y_im[I, :] .= (w[I, :, i] .* s.Cuf[I, :, i] .* Ti) .* (1 .- S[I, :]) .+
                                 (w[I, :, i] .* s.Cu[I, :, i] .* Ti) .* S[I, :]
            lC = view(l.Ct, :, :, i)
            @views y_ex[I, :] .= (l.w[I, :, i] .* l.Cuf[I, :, i] .* Ti) .* (1 .- S[I, :]) .+
                                 (l.w[I, :, i] .* l.Cu[I, :, i] .* Ti) .* S[I, :] .-
                                 (sgs[I, :] .* Cs[I, :] .+ sgn[I, :] .* Cn[I, :] .+ Ti) .*
                                 lC[I, :] .+
                                 ixs[I, :] .* Cs_plus[I, :] .* lC[1:(end - 2), :] .-
                                 (1 .- ixs[I, :]) .* Cs_min[I, :] .* lC[3:end, :] .+
                                 ixn[I, :] .* Cn_plus[I, :] .* circshift(lC[I, :], (0, 1)) .-
                                 (1 .- ixn[I, :]) .* Cn_min[I, :] .*
                                 circshift(lC[I, :], (0, -1))
            y_i = zeros(nxp, nyp)
            @views y_i[I, :] .= lC[I, :] .+ alpha .* y_im[I, :] .+ (1 - alpha) .* y_ex[I, :]
            bo == :flux && (y_i[1, :] .= p.offshore_flux .* s.Cu0[1, :, i])
            bn == :flux && (y_i[end, :] .= p.onshore_flux .* s.Cu0[end, :, i])
            bo == :constant && (y_i[1, :] .= p.constant_offshore_flux ./ s.u[1, :, i])
            bn == :constant && (y_i[end, :] .= p.constant_onshore_flux ./ s.u[end, :, i])
            rhs = vec(y_i)
            if cache === nothing
                cache = _linsolve_cache(A, rhs, linsolve)
            end
            Ct_i = _linsolve!(cache, rhs)
            prevent_tiny_negatives!(Ct_i, p.max_error)
            if minimum(Ct_i) < 0
                neg = Ct_i .< 0
                pos = .!neg
                spos = sum(Ct_i[pos])
                if spos != 0
                    Ct_i[pos] .*= 1 + sum(Ct_i[neg]) / spos
                else
                    Ct_i[pos] .= 0
                end
                Ct_i[neg] .= 0
            end
            w_i = vec(w[:, :, i])
            pickup_i = (w_i .* Cu_i .- Ct_i) ./ p.T .* dt
            deficit = pickup_i .- mass_i
            ixc = (deficit .> p.max_error) .& (w_i .* Cu_i .> 0)
            if !any(ixc)
                pickup_i = min.(pickup_i, mass_i)
                conv = true
                break
            else
                w_i[ixc] .= (mass_i[ixc] .* p.T ./ dt .+ Ct_i[ixc]) ./ Cu_i[ixc]
                w[:, :, i] .= reshape(w_i, nxp, nyp)
            end
        end
        conv ||
            @warn "iteration not converged" nrcells = count(ixc) fraction = i time = m.time
        Ct[:, :, i] .= reshape(Ct_i, nxp, nyp)
        pickup[:, :, i] .= reshape(pickup_i, nxp, nyp)
    end
    any(1 .- sum(w; dims=3) .> p.max_error) && _count!(m, :supplylim)
    return (; Ct, qs=Ct .* s.us, qn=Ct .* s.un, pickup, w, w_init, w_air, w_bed)
end

"""
    solve_euler_forward(s, m, dt)

The explicit AeoLiS solver (`solve_EF`): first-order upwind advection on a grid with ghost
cells (the first cross-shore column and the ghost column are empty, the lateral edges wrap)
plus the pickup limited by the available sediment.
"""
function solve_euler_forward(s::AeolisState, m::AeolisModel, dt::Real)
    p = m.p
    nxp, nyp = m.nx + 1, m.ny + 1
    Ct = copy(s.Ct)
    pickup = copy(s.pickup)
    Ts = p.T
    w_init, w_air, w_bed = compute_weights(s, m)
    w = m.time == 0 ? _initial_weights(s, m, w_init) : copy(w_init)
    us0 = s.us[:, :, 1]
    un0 = s.un[:, :, 1]
    ix1 = Float64.((us0 .>= 0) .& (un0 .>= 0))
    ix2 = Float64.((us0 .< 0) .& (un0 .>= 0))
    ix3 = Float64.((us0 .< 0) .& (un0 .< 0))
    ix4 = Float64.((us0 .> 0) .& (un0 .< 0))
    for i in 1:(m.nf)
        Cg = zeros(nxp + 2, nyp + 2)
        Cg[2:(end - 1), 2:(end - 1)] .= Ct[:, :, i]
        Cg[1:2, :] .= 0
        Cg[:, 1] .= Cg[:, end - 1]
        Cg[:, end] .= Cg[:, 2]
        @views pickup[:, :, i] .= dt .*
                                  (min.(s.Cu[:, :, i], s.mass[:, :, 1, i] .+ Ct[:, :, i]) .-
                                   Ct[:, :, i]) ./ Ts
        C = Cg[2:(end - 1), 2:(end - 1)]
        W = Cg[1:(end - 2), 2:(end - 1)]
        E = Cg[3:end, 2:(end - 1)]
        S = Cg[2:(end - 1), 1:(end - 2)]
        N = Cg[2:(end - 1), 3:end]
        us = s.us[:, :, i]
        un = s.un[:, :, i]
        ds, dn = s.ds, s.dn
        new = @. C + ix1 * (-dt * us * (C - W) / ds - dt * un * (C - S) / dn) +
                 ix2 * (+dt * us * (C - E) / ds - dt * un * (C - S) / dn) +
                 ix3 * (+dt * us * (C - E) / ds + dt * un * (C - N) / dn) +
                 ix4 * (-dt * us * (C - W) / ds + dt * un * (C - N) / dn) + pickup[:, :, i]
        Ct[:, :, i] .= new
    end
    return (; Ct, qs=Ct .* s.us, qn=Ct .* s.un, pickup, w, w_init, w_air, w_bed)
end

"""
    solve_pieter(s, m, dt, alpha, beta, linsolve)

The finite-volume AeoLiS solver (`solve_pieter`): upwind fluxes on cell faces, the change of
concentration solved for iteratively with erosion limited by the bed, and flux, constant,
gradient or circular boundaries. At the first step (`t == 0`) only the weights are set.
"""
function solve_pieter(s::AeolisState, m::AeolisModel, dt::Real, alpha::Real, beta::Real,
                      linsolve)
    p = m.p
    l = m.prev
    nxp, nyp = m.nx + 1, m.ny + 1
    nf = m.nf
    Ct = copy(s.Ct)
    pickup = copy(s.pickup)
    Ts = p.T
    w_init, w_air, w_bed = compute_weights(s, m)
    if m.time == 0
        return (; w=_initial_weights(s, m, w_init))
    end
    w = copy(w_init)
    bo, bn, bl = p.boundary_offshore, p.boundary_onshore, p.boundary_lateral
    for i in 1:nf
        us = s.us[:, :, i]
        un = s.un[:, :, i]
        ufs = zeros(nxp + 1, nyp)
        ufn = zeros(nxp, nyp + 1)
        ufs[2:(end - 1), :] .= 0.5 .* us[1:(end - 1), :] .+ 0.5 .* us[2:end, :]
        ufn[:, 2:(end - 1)] .= 0.5 .* un[:, 1:(end - 1)] .+ 0.5 .* un[:, 2:end]
        ufs[1, :] .= us[1, :]
        ufs[end, :] .= us[end, :]
        if bl == :circular
            ufn[:, 1] .= 0.5 .* un[:, 1] .+ 0.5 .* un[:, end]
            ufn[:, end] .= ufn[:, 1]
        else
            ufn[:, 1] .= un[:, 1]
            ufn[:, end] .= un[:, end]
        end
        b = abs(beta)
        ixfs = b >= 1 ? Float64.(ufs .>= 0) : fill(b, size(ufs))
        ixfn = b >= 1 ? Float64.(ufn .>= 0) : fill(b, size(ufn))
        dsdn, ds, dn = s.dsdn, s.ds, s.dn
        A0 = zeros(nxp, nyp)
        Apx = zeros(nxp, nyp)
        Ap1 = zeros(nxp, nyp)
        Amx = zeros(nxp, nyp)
        Am1 = zeros(nxp, nyp)
        A0 .+= dsdn ./ dt
        A0 .+= dsdn ./ Ts .* alpha
        F = ufs[2:(end - 1), :]
        G = ixfs[2:(end - 1), :]
        @views A0[2:end, :] .-= dn[2:end, :] .* F .* (1 .- G) .* alpha
        @views Am1[2:end, :] .-= dn[2:end, :] .* F .* G .* alpha
        @views A0[1:(end - 1), :] .+= dn[1:(end - 1), :] .* F .* G .* alpha
        @views Ap1[1:(end - 1), :] .+= dn[1:(end - 1), :] .* F .* (1 .- G) .* alpha
        Fn = ufn[:, 2:(end - 1)]
        Gn = ixfn[:, 2:(end - 1)]
        @views A0[:, 2:end] .-= ds[:, 2:end] .* Fn .* (1 .- Gn) .* alpha
        @views Amx[:, 2:end] .-= ds[:, 2:end] .* Fn .* Gn .* alpha
        @views A0[:, 1:(end - 1)] .+= ds[:, 1:(end - 1)] .* Fn .* Gn .* alpha
        @views Apx[:, 1:(end - 1)] .+= ds[:, 1:(end - 1)] .* Fn .* (1 .- Gn) .* alpha
        if bo == :constant
            A0[1, :] .= 1
            Apx[1, :] .= 0
            Amx[1, :] .= 0
            Ap1[1, :] .= 0
            Am1[1, :] .= 0
        elseif bo == :gradient
            @views A0[1, :] .-= dn[1, :] .* ufs[2, :] .* ixfs[2, :] .* alpha
            @views Ap1[1, :] .-= dn[1, :] .* ufs[2, :] .* (1 .- ixfs[2, :]) .* alpha
        elseif bo != :flux
            throw(UnsupportedOptionError(:boundary_offshore, bo,
                                         "the pieter solver supports flux, constant and gradient offshore boundaries"))
        end
        if bn == :constant
            A0[end, :] .= 1
            Apx[end, :] .= 0
            Amx[end, :] .= 0
            Ap1[end, :] .= 0
            Am1[end, :] .= 0
        elseif bn == :gradient
            @views A0[end, :] .+= dn[end, :] .* ufs[end - 1, :] .* (1 .- ixfs[end - 1, :]) .*
                                  alpha
            @views Am1[end, :] .+= dn[end, :] .* ufs[end - 1, :] .* ixfs[end - 1, :] .*
                                   alpha
        elseif bn != :flux
            throw(UnsupportedOptionError(:boundary_onshore, bn,
                                         "the pieter solver supports flux, constant and gradient onshore boundaries"))
        end
        if bl == :constant
            for c in (1, nyp)
                A0[:, c] .= 1
                Apx[:, c] .= 0
                Amx[:, c] .= 0
                Ap1[:, c] .= 0
                Am1[:, c] .= 0
            end
        elseif bl == :gradient
            @views A0[:, 1] .-= ds[:, 1] .* ufn[:, 2] .* ixfn[:, 2] .* alpha
            @views Apx[:, 1] .-= ds[:, 1] .* ufn[:, 2] .* (1 .- ixfn[:, 2]) .* alpha
            @views A0[:, end] .+= ds[:, end] .* ufn[:, end - 1] .* (1 .- ixfn[:, end - 1]) .*
                                  alpha
            @views Amx[:, end] .+= ds[:, end] .* ufn[:, end - 1] .* ixfn[:, end - 1] .*
                                   alpha
        elseif bl == :circular
            @views A0[:, 1] .-= ds[:, 1] .* ufn[:, 1] .* (1 .- ixfn[:, 1]) .* alpha
            @views Amx[:, 1] .-= ds[:, 1] .* ufn[:, 1] .* ixfn[:, 1] .* alpha
            @views A0[:, end] .+= ds[:, end] .* ufn[:, end] .* ixfn[:, end] .* alpha
            @views Apx[:, end] .+= ds[:, end] .* ufn[:, end] .* (1 .- ixfn[:, end]) .* alpha
        elseif bl != :flux
            throw(UnsupportedOptionError(:boundary_lateral, bl,
                                         "the pieter solver supports flux, constant, gradient and circular lateral boundaries"))
        end
        n = nxp * nyp
        j = nxp
        ny = m.ny
        A = diags((vec(Apx)[1:j], vec(Amx)[(j + 1):end], vec(Am1)[2:end], vec(A0),
                   vec(Ap1)[1:(end - 1)], vec(Apx)[(j + 1):end], vec(Amx)[1:j]),
                  (-j * ny, -j, -1, 0, 1, j, j * ny), n)
        renormalize_weights!(w, i)
        cache = nothing
        conv = false
        ixc = falses(n)
        mass_i = vec(s.mass[:, :, 1, i])
        Cu_i = vec(s.Cu[:, :, i])
        lC = l.Ct[:, :, i]
        for it in 1:(p.max_iter)
            _count!(m, :matrixsolve)
            C = Ct[:, :, i]
            Cxs = zeros(nxp + 1, nyp)
            Cxn = zeros(nxp, nyp + 1)
            @views Cxs[2:(end - 1), :] .= G .* (alpha .* C[1:(end - 1), :] .+
                                                (1 - alpha) .* lC[1:(end - 1), :]) .+
                                          (1 .- G) .* (alpha .* C[2:end, :] .+
                                                       (1 - alpha) .* lC[2:end, :])
            @views Cxn[:, 2:(end - 1)] .= Gn .* (alpha .* C[:, 1:(end - 1)] .+
                                                 (1 - alpha) .* lC[:, 1:(end - 1)]) .+
                                          (1 .- Gn) .* (alpha .* C[:, 2:end] .+
                                                        (1 - alpha) .* lC[:, 2:end])
            if bl == :circular
                @views Cxn[:, 1] .= ixfn[:, 1] .*
                                    (alpha .* C[:, end] .+ (1 - alpha) .* lC[:, end]) .+
                                    (1 .- ixfn[:, 1]) .*
                                    (alpha .* C[:, 1] .+ (1 - alpha) .* lC[:, 1])
                Cxn[:, end] .= Cxn[:, 1]
            end
            D_i = dsdn ./ Ts .* (alpha .* C .+ (1 - alpha) .* lC)
            A_i = dsdn ./ Ts .* s.mass[:, :, 1, i] .+ D_i
            U_i = dsdn ./ Ts .* (w[:, :, i] .* alpha .* s.Cu[:, :, i] .+
                                 (1 - alpha) .* l.w[:, :, i] .* l.Cu[:, :, i])
            E_i = min.(U_i, A_i)
            y = zeros(nxp, nyp)
            y .-= dsdn ./ dt .* (C .- lC)
            y .+= E_i .- D_i
            @views y[2:end, :] .+= dn[2:end, :] .* F .* Cxs[2:(end - 1), :]
            @views y[1:(end - 1), :] .-= dn[1:(end - 1), :] .* F .* Cxs[2:(end - 1), :]
            @views y[:, 2:end] .+= ds[:, 2:end] .* Fn .* Cxn[:, 2:(end - 1)]
            @views y[:, 1:(end - 1)] .-= ds[:, 1:(end - 1)] .* Fn .* Cxn[:, 2:(end - 1)]
            if bo == :flux
                @views y[1, :] .+= dn[1, :] .* ufs[1, :] .* s.Cu0[1, :, i] .*
                                   p.offshore_flux
            elseif bo == :constant
                y[1, :] .= 0
            elseif bo == :gradient
                @views y[1, :] .+= dn[2, :] .* ufs[2, :] .* Cxs[2, :]
            end
            if bn == :flux
                @views y[end, :] .+= dn[end, :] .* ufs[end, :] .* s.Cu0[end, :, i] .*
                                     p.onshore_flux
            elseif bn == :constant
                y[end, :] .= 0
            elseif bn == :gradient
                @views y[end, :] .-= dn[end - 1, :] .* ufs[end - 1, :] .* Cxs[end - 1, :]
            end
            if bl == :flux
                @views y[:, 1] .+= ds[:, 1] .* ufn[:, 1] .* s.Cu0[:, 1, i] .* p.lateral_flux
                @views y[:, end] .-= ds[:, end] .* ufn[:, end] .* s.Cu0[:, end, i] .*
                                     p.lateral_flux
            elseif bl == :constant
                y[:, 1] .= 0
                y[:, end] .= 0
            elseif bl == :gradient
                @views y[:, end] .-= ds[:, end - 1] .* ufn[:, end - 1] .* Cxn[:, end - 1]
                @views y[:, 1] .+= ds[:, 2] .* ufn[:, 2] .* Cxn[:, 2]
            elseif bl == :circular
                @views y[:, 1] .+= ds[:, 1] .* ufn[:, 1] .* Cxn[:, 1]
                @views y[:, end] .-= ds[:, end] .* ufn[:, end] .* Cxn[:, end]
            end
            rhs = vec(y)
            if cache === nothing
                cache = _linsolve_cache(A, rhs, linsolve)
            end
            Ct_i = vec(C) .+ _linsolve!(cache, rhs)
            prevent_tiny_negatives!(Ct_i, p.max_error)
            Ct_i[Ct_i .< 0] .= 0
            w_i = vec(w[:, :, i])
            pickup_i = (w_i .* Cu_i .- Ct_i) ./ Ts .* dt
            deficit = pickup_i .- mass_i
            ixc = (deficit .> p.max_error) .& (w_i .* Cu_i .> 0)
            # AeoLiS stores the pickup before limiting it by the bed
            pickup[:, :, i] .= reshape(pickup_i, nxp, nyp)
            Ct[:, :, i] .= reshape(Ct_i, nxp, nyp)
            if !any(ixc)
                conv = true
                break
            else
                w_i[ixc] .= (mass_i[ixc] .* Ts ./ dt .+ Ct_i[ixc]) ./ Cu_i[ixc]
                w[:, :, i] .= reshape(w_i, nxp, nyp)
            end
        end
        conv ||
            @warn "iteration not converged" nrcells = count(ixc) fraction = i time = m.time
    end
    any(1 .- sum(w; dims=3) .> p.max_error) && _count!(m, :supplylim)
    qs = Ct .* s.us
    qn = Ct .* s.un
    return (; Ct, qs, qn, pickup, w, w_init, w_air, w_bed, q=hypot.(qs, qn))
end

"""
    solve_steady_state(s, m, dt)

The steady-state AeoLiS solver (`solve_SS` with `sweep`): the concentration in equilibrium
with the current wind, found by sweeping the grid in the four upwind directions, with the
pickup limited by the top layer; circular boundaries are iterated to convergence.
"""
function solve_steady_state(s::AeolisState, m::AeolisModel, dt::Real)
    p = m.p
    w_init, w_air, w_bed = compute_weights(s, m)
    w = m.time == 0 ? _initial_weights(s, m, w_init) : copy(w_init)
    Ct = zeros(size(s.Ct))
    p.boundary_offshore == :flux && (Ct[1, :, 1] .= s.Cu0[1, :, 1])
    p.boundary_onshore == :flux && (Ct[end, :, 1] .= s.Cu0[end, :, 1])
    if p.boundary_offshore == :circular
        Ct[1, :, 1] .= -1
        Ct[end, :, 1] .= -1
    end
    if p.boundary_offshore == :re_circular
        Ct[1, :, 1] .= -2
        Ct[end, :, 1] .= -2
    end
    if p.boundary_lateral == :circular
        Ct[:, 1, 1] .= -1
        Ct[:, end, 1] .= -1
    end
    if p.boundary_lateral == :re_circular
        Ct[:, 1, 1] .= -2
        Ct[:, end, 1] .= -2
    end
    Ct, pickup = sweep!(Ct, copy(s.Cu), s.mass, dt, p.T, s.ds, s.dn, s.us, s.un, w)
    qs = Ct .* s.us
    qn = Ct .* s.un
    return (; Ct, qs, qn, pickup, w, w_init, w_air, w_bed, q=hypot.(qs, qn))
end

"""
    sweep!(Ct, Cu, mass, dt, Ts, ds, dn, us, un, w)

AeoLiS's steady-state sweep. Boundary types are passed in `Ct`'s first fraction as in
AeoLiS (`-1` on the edges for circular, `-2` for re-circulating). Indices below follow the
Python code (0-based `n` along y, `s` along x) and are converted at each access.
"""
function sweep!(Ct, Cu, mass, dt, Ts, ds, dn, us, un, w)
    nxp, nyp, nf = size(Cu)
    pickup = zeros(size(Cu))
    circ_lateral = Ct[2, 1, 1] == -1
    if circ_lateral
        Ct[:, 1, 1] .= 0
        Ct[:, end, 1] .= 0
    end
    circ_offshore = Ct[1, 2, 1] == -1
    if circ_offshore
        Ct[1, :, 1] .= 0
        Ct[end, :, 1] .= 0
    end
    recirc_offshore = Ct[1, 2, 1] == -2
    if recirc_offshore
        Ct[1, :, 1] .= 0
        Ct[end, :, 1] .= 0
    end
    ufs = zeros(nxp + 1, nyp, nf)
    ufn = zeros(nxp, nyp + 1, nf)
    ufs[2:(end - 1), :, :] .= 0.5 .* us[1:(end - 1), :, :] .+ 0.5 .* us[2:end, :, :]
    ufn[:, 2:(end - 1), :] .= 0.5 .* un[:, 1:(end - 1), :] .+ 0.5 .* un[:, 2:end, :]
    ufs[1, :, :] .= ufs[2, :, :]
    ufs[end, :, :] .= ufs[end - 1, :, :]
    ufn[:, 1, :] .= ufn[:, 2, :]
    ufn[:, end, :] .= ufn[:, end - 1, :]
    ufs[1, :, :] .= (ufs[1, :, :] .+ ufs[end, :, :]) ./ 2
    ufs[end, :, :] .= ufs[1, :, :]
    ufs[:, 1, :] .= (ufs[:, 1, :] .+ ufs[:, end, :]) ./ 2
    ufs[:, end, :] .= ufs[:, 1, :]
    ufn[1, :, :] .= (ufn[1, :, :] .+ ufn[end, :, :]) ./ 2
    ufn[end, :, :] .= ufn[1, :, :]
    ufn[:, 1, :] .= (ufn[:, 1, :] .+ ufn[:, end, :]) ./ 2
    ufn[:, end, :] .= ufn[:, 1, :]
    k = 0
    Ct_last = copy(Ct)
    while k == 0 || any(abs.(Ct[:, :, 1] .- Ct_last[:, :, 1]) .> 1e-10)
        Ct_last = copy(Ct)
        if circ_lateral
            a, b = copy(Ct[:, end, 1]), copy(Ct[:, 1, 1])
            Ct[:, 1, 1] .= a
            Ct[:, end, 1] .= b
        end
        if circ_offshore
            a, b = copy(Ct[end, :, 1]), copy(Ct[1, :, 1])
            Ct[1, :, 1] .= a
            Ct[end, :, 1] .= b
        end
        if recirc_offshore
            a, b = mean(Ct[end - 1, :, 1]), mean(Ct[2, :, 1])
            Ct[1, :, 1] .= a
            Ct[end, :, 1] .= b
        end
        visited = falses(nxp, nyp)
        _sweep_quadrants!(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, nf)
        pickup[:, 1, 1] .= pickup[:, 2, 1]
        pickup[:, end, 1] .= pickup[:, end - 1, 1]
        k += 1
    end
    return Ct, pickup
end

function _sweep_quadrants!(Ct, Cu, mass, pickup, dt, Ts, ds, dn, ufs, ufn, w, visited, nf)
    nxp, nyp = size(Cu, 1), size(Cu, 2)
    # 0-based accessors with AeoLiS's (n, s) order
    C(n, s, f) = Ct[s + 1, n + 1, f + 1]
    setC!(v, n, s, f) = (Ct[s + 1, n + 1, f + 1] = v)
    FS(n, s, f) = ufs[s + 1, n + 1, f + 1]
    FN(n, s, f) = ufn[s + 1, n + 1, f + 1]
    DS(n, s) = ds[s + 1, n + 1]
    DN(n, s) = dn[s + 1, n + 1]
    WC(n, s, f) = w[s + 1, n + 1, f + 1] * Cu[s + 1, n + 1, f + 1]
    M(n, s, f) = mass[s + 1, n + 1, 1, f + 1]
    V(n, s) = visited[s + 1, n + 1]
    setV!(n, s) = (visited[s + 1, n + 1] = true)
    P(n, s, f) = pickup[s + 1, n + 1, f + 1]
    setP!(v, n, s, f) = (pickup[s + 1, n + 1, f + 1] = v)
    NY, NX = nyp, nxp
    function finish!(n, s, f, num_lim_terms, den_lim)
        return nothing
    end
    # quadrant 1: wind towards +x and +y
    for n in 1:(NY - 1), s in 1:(NX - 1)
        if !V(n, s) && FN(n, s, 0) >= 0 && FS(n, s, 0) >= 0 && FN(n + 1, s, 0) >= 0 &&
           FS(n, s + 1, 0) >= 0
            for f in 0:(nf - 1)
                num = C(n - 1, s, f) * FN(n, s, f) * DS(n, s) +
                      C(n, s - 1, f) * FS(n, s, f) * DN(n, s) +
                      WC(n, s, f) * DS(n, s) * DN(n, s) / Ts
                den = FN(n + 1, s, f) * DS(n, s) + FS(n, s + 1, f) * DN(n, s) +
                      DS(n, s) * DN(n, s) / Ts
                setC!(num / den, n, s, f)
                setP!((WC(n, s, f) - C(n, s, f)) * dt / Ts, n, s, f)
                if P(n, s, f) > M(n, s, f)
                    setP!(M(n, s, f), n, s, f)
                    num = C(n - 1, s, f) * FN(n, s, f) * DS(n, s) +
                          C(n, s - 1, f) * FS(n, s, f) * DN(n, s) +
                          P(n, s, f) * DS(n, s) * DN(n, s) / dt
                    den = FN(n + 1, s, f) * DS(n, s) + FS(n, s + 1, f) * DN(n, s)
                    setC!(num / den, n, s, f)
                end
            end
            setV!(n, s)
        end
    end
    # quadrant 2: towards -x and +y
    for n in 1:(NY - 1), s in (NX - 2):-1:0
        if !V(n, s) && FN(n, s, 0) >= 0 && FS(n, s, 0) <= 0 && FN(n + 1, s, 0) >= 0 &&
           FS(n, s + 1, 0) <= 0
            for f in 0:(nf - 1)
                num = C(n - 1, s, f) * FN(n, s, f) * DS(n, s) +
                      -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s) +
                      WC(n, s, f) * DS(n, s) * DN(n, s) / Ts
                den = FN(n + 1, s, f) * DS(n, s) + -FS(n, s, f) * DN(n, s) +
                      DS(n, s) * DN(n, s) / Ts
                setC!(num / den, n, s, f)
                setP!((WC(n, s, f) - C(n, s, f)) * dt / Ts, n, s, f)
                if P(n, s, f) > M(n, s, f)
                    setP!(M(n, s, f), n, s, f)
                    num = C(n - 1, s, f) * FN(n, s, f) * DS(n, s) +
                          -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s) +
                          P(n, s, f) * DS(n, s) * DN(n, s) / dt
                    den = FN(n + 1, s, f) * DS(n, s) + -FS(n, s, f) * DN(n, s)
                    setC!(num / den, n, s, f)
                end
            end
            setV!(n, s)
        end
    end
    # quadrant 3: towards -x and -y (AeoLiS uses dn for the y-faces here)
    for n in (NY - 2):-1:0, s in (NX - 2):-1:0
        if !V(n, s) && FN(n, s, 0) <= 0 && FS(n, s, 0) <= 0 && FN(n + 1, s, 0) <= 0 &&
           FS(n, s + 1, 0) <= 0
            for f in 0:(nf - 1)
                num = -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s) +
                      -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s) +
                      WC(n, s, f) * DS(n, s) * DN(n, s) / Ts
                den = -FN(n, s, f) * DN(n, s) + -FS(n, s, f) * DN(n, s) +
                      DS(n, s) * DN(n, s) / Ts
                setC!(num / den, n, s, f)
                setP!((WC(n, s, f) - C(n, s, f)) * dt / Ts, n, s, f)
                if P(n, s, f) > M(n, s, f)
                    setP!(M(n, s, f), n, s, f)
                    num = -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s) +
                          -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s) +
                          P(n, s, f) * DS(n, s) * DN(n, s) / dt
                    den = -FN(n, s, f) * DN(n, s) + -FS(n, s, f) * DN(n, s)
                    setC!(num / den, n, s, f)
                end
            end
            setV!(n, s)
        end
    end
    # quadrant 4: towards +x and -y (AeoLiS uses dn for the y-faces here)
    for n in (NY - 2):-1:0, s in 1:(NX - 1)
        if !V(n, s) && FN(n, s, 0) <= 0 && FS(n, s, 0) >= 0 && FN(n + 1, s, 0) <= 0 &&
           FS(n, s + 1, 0) >= 0
            for f in 0:(nf - 1)
                num = C(n, s - 1, f) * FS(n, s, f) * DN(n, s) +
                      -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s) +
                      WC(n, s, f) * DS(n, s) * DN(n, s) / Ts
                den = FS(n, s + 1, f) * DN(n, s) + -FN(n, s, f) * DN(n, s) +
                      DS(n, s) * DN(n, s) / Ts
                setC!(num / den, n, s, f)
                setP!((WC(n, s, f) - C(n, s, f)) * dt / Ts, n, s, f)
                if P(n, s, f) > M(n, s, f)
                    setP!(M(n, s, f), n, s, f)
                    num = C(n, s - 1, f) * FS(n, s, f) * DN(n, s) +
                          -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s) +
                          P(n, s, f) * DS(n, s) * DN(n, s) / dt
                    den = FS(n, s + 1, f) * DN(n, s) + -FN(n, s, f) * DN(n, s)
                    setC!(num / den, n, s, f)
                end
            end
            setV!(n, s)
        end
    end
    # the remaining cells, with a generic upwind stencil
    for n in (NY - 2):-1:0, s in 1:(NX - 1)
        if !V(n, s) && n != 0 && s != NX - 1
            for f in 0:(nf - 1)
                num = WC(n, s, f) * DS(n, s) * DN(n, s) / Ts
                den = DS(n, s) * DN(n, s) / Ts
                FN(n, s, 0) > 0 && (num += C(n - 1, s, f) * FN(n, s, f) * DS(n, s))
                FS(n, s, 0) > 0 && (num += C(n, s - 1, f) * FS(n, s, f) * DN(n, s))
                if FN(n + 1, s, 0) < 0
                    num += -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s)
                elseif FN(n + 1, s, 0) > 0
                    den += FN(n + 1, s, f) * DS(n, s)
                end
                if FS(n, s + 1, 0) < 0
                    num += -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s)
                elseif FS(n, s + 1, 0) > 0
                    den += FS(n, s + 1, f) * DN(n, s)
                end
                FN(n, s, 0) < 0 && (den += -FN(n, s, f) * DN(n, s))
                FS(n, s, 0) < 0 && (den += -FS(n, s, f) * DN(n, s))
                setC!(num / den, n, s, f)
                setP!((WC(n, s, f) - C(n, s, f)) * dt / Ts, n, s, f)
                if P(n, s, f) > M(n, s, f)
                    setP!(M(n, s, f), n, s, f)
                    nl = P(n, s, f) * DS(n, s) * DN(n, s) / dt
                    dl = 0.0
                    FN(n, s, 0) > 0 && (nl += C(n - 1, s, f) * FN(n, s, f) * DS(n, s))
                    FS(n, s, 0) > 0 && (nl += C(n, s - 1, f) * FS(n, s, f) * DN(n, s))
                    if FN(n + 1, s, 0) < 0
                        nl += -C(n + 1, s, f) * FN(n + 1, s, f) * DN(n, s)
                    elseif FN(n + 1, s, 0) > 0
                        dl += FN(n + 1, s, f) * DS(n, s)
                    end
                    if FS(n, s + 1, 0) < 0
                        nl += -C(n, s + 1, f) * FS(n, s + 1, f) * DN(n, s)
                    elseif FS(n, s + 1, 0) > 0
                        dl += FS(n, s + 1, f) * DN(n, s)
                    end
                    FN(n, s, 0) < 0 && (dl += -FN(n, s, f) * DN(n, s))
                    FS(n, s, 0) < 0 && (dl += -FS(n, s, f) * DN(n, s))
                    setC!(nl / dl, n, s, f)
                end
            end
            setV!(n, s)
        end
    end
    return Ct
end
