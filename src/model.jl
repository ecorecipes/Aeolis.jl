# One time step of AeoLiS (`AeoLiS.update`): wind and shear, fences and vegetation, the time
# step, hydrodynamics, threshold, equilibrium transport, the transport solver, the bed update,
# avalanching, bed reset, dune erosion and vegetation growth.

"""
    set_timestep!(m, s, dt, scheme)

The time step of the next update: `dt` when positive, otherwise the `dt` parameter, limited by
the CFL condition for the explicit scheme and by `max_bedlevel_change`. Returns `false` for
`dt < -1`, which skips `|dt|` seconds without computing.
"""
function set_timestep!(m::AeolisModel, s::AeolisState, dt::Real, scheme::Symbol)
    p = m.p
    if dt > 0
        step = Float64(dt)
    elseif dt < -1
        return false, Float64(dt)
    else
        step = p.dt
    end
    if scheme == :euler_forward && p.CFL > 0
        dtref = maximum(abs.(s.uws) ./ s.ds) + maximum(abs.(s.uwn) ./ s.dn)
        dtref > 0 && (step = min(step, p.CFL / dtref))
    end
    mz = maximum(s.dzb)
    if p.max_bedlevel_change != 999.0 && mz != 0 && m.dt_prev != 0
        step = min(step, m.dt_prev * p.max_bedlevel_change / mz)
    end
    m.dt_opt = step
    return true, step
end

function _apply_transport!(s::AeolisState, r::NamedTuple)
    for (k, v) in pairs(r)
        getfield(s, k) .= v
    end
    return s
end

function _transport!(s, m, dt, scheme::Symbol, solver::Symbol, linsolve)
    alpha = scheme == :euler_forward ? 0.0 : scheme == :euler_backward ? 1.0 : 0.5
    r = if solver == :trunk
        scheme == :euler_forward ? solve_euler_forward(s, m, dt) :
        solve_trunk(s, m, dt, alpha, 1.0, linsolve)
    elseif solver == :pieter
        solve_pieter(s, m, dt, alpha, 1.0, linsolve)
    elseif solver == :steadystate
        # AeoLiS: the sweep with euler_forward and euler_backward, solve_steadystate with
        # crank_nicolson
        scheme == :crank_nicolson ? solve_trunk(s, m, dt, 1.0, 1.0, linsolve; steady=true) :
        solve_steady_state(s, m, dt)
    elseif solver == :steadystatepieter
        solve_pieter(s, m, dt, 1.0, 1.0, linsolve; steady=true)
    else
        throw(UnsupportedOptionError(:solver, solver,
                                     "use trunk, pieter, steadystate or steadystatepieter"))
    end
    return _apply_transport!(s, r)
end

# a transport solver passed as an object (SciMLSteadyState)
function _transport!(s, m, dt, scheme::Symbol, solver, linsolve)
    return _apply_transport!(s, solve_sciml_steady_state(s, m, dt, solver))
end

"""
    aeolis_step!(s, m, t, dt, scheme, solver, linsolve)

Advance the state `s` from time `t` by one AeoLiS update with requested step `dt` (`-1` for
the configured step). `solver` is an AeoLiS solver name (`:trunk`, `:pieter`, `:steadystate`,
`:steadystatepieter`) or a [`SciMLSteadyState`](@ref). Returns the time step used, or
`nothing` when no update was made.
"""
function aeolis_step!(s::AeolisState, m::AeolisModel, t::Real, dt::Real, scheme::Symbol,
                      solver, linsolve)
    p = m.p
    m.time = t
    copyto_state!(m.prev, s)
    interpolate_wind!(s, m, t)
    grid_rotate!(s, m, m.alpha)
    sum(s.uw) != 0 && compute_shear!(s, m)
    p.process_fences && update_fences!(s, m)
    if p.process_vegetation
        p.veg_model == :marram && marram_cover!(s, m)
        vegetation_shear!(s, m)
    end
    m.dt_prev = m.dt_opt
    ok, step = set_timestep!(m, s, dt, scheme)
    if !ok
        grid_rotate!(s, m, -m.alpha)
        return nothing
    end
    interpolate_hydro!(s, m, t)
    update_hydro!(s, m, step, t)
    mixtoplayer!(s, m)
    compute_threshold!(s, m)
    equilibrium!(s, m)
    _transport!(s, m, step, scheme, solver, linsolve)
    update_bed!(s, m)
    angle_of_repose!(s, m)
    avalanche!(s, m)
    wet_bed_reset!(s, m)
    average_change!(s, m)
    if p.process_dune_erosion
        dune_erosion!(s, m, t)
        angle_of_repose!(s, m)
        avalanche!(s, m)
    end
    if p.process_vegetation
        if p.veg_model == :marram
            marram_update!(s, m, t)
        else
            germinate!(s, m)
            grow_vegetation!(s, m)
        end
    end
    _count!(m, :time)
    grid_rotate!(s, m, -m.alpha)
    return step
end
