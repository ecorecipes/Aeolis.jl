# The SciML interface: `AeolisProblem`, the numerical schemes as algorithms, an integrator
# with `init`, `step!`, `solve!`, `reinit!` and SciMLBase callbacks, and `AeolisSolution`.

"""
    AbstractTransportSolver

Supertype of the transport solvers used inside a scheme: [`Trunk`](@ref),
[`Pieter`](@ref) and [`SteadyState`](@ref).
"""
abstract type AbstractTransportSolver end

"""
    Trunk()

The default AeoLiS transport solver: an upwind discretisation solved as a sparse linear
system per fraction (`solver = trunk`).
"""
struct Trunk <: AbstractTransportSolver end

"""
    Pieter()

The finite-volume AeoLiS transport solver with fluxes on cell faces (`solver = pieter`).
"""
struct Pieter <: AbstractTransportSolver end

"""
    SteadyState()

The steady-state AeoLiS transport solver, sweeping the grid upwind (`solver = steadystate`).
"""
struct SteadyState <: AbstractTransportSolver end

_solver_symbol(::Trunk) = :trunk
_solver_symbol(::Pieter) = :pieter
_solver_symbol(::SteadyState) = :steadystate
function _solver_from_symbol(s::Symbol)
    s == :trunk && return Trunk()
    s == :pieter && return Pieter()
    s == :steadystate && return SteadyState()
    return throw(UnsupportedOptionError(:solver, s, "the $(s) solver is not ported"))
end

"""
    AbstractAeolisAlgorithm

Supertype of the numerical schemes: [`EulerBackward`](@ref), [`CrankNicolson`](@ref) and
[`EulerForward`](@ref), each with a transport solver and an optional LinearSolve.jl
algorithm for the sparse systems.
"""
abstract type AbstractAeolisAlgorithm <: SciMLBase.AbstractSciMLAlgorithm end

for (name, sym, doc) in
    ((:EulerBackward, :euler_backward, "implicit Euler (AeoLiS default)"),
     (:CrankNicolson, :crank_nicolson, "Crank-Nicolson"),
     (:EulerForward, :euler_forward,
      "explicit Euler, with the time step limited by the CFL number"))
    @eval begin
        """
            $($(string(name)))(solver=Trunk(); linsolve=nothing)

        The $($doc) scheme with transport solver `solver`. `linsolve` is a LinearSolve.jl
        algorithm for the sparse systems (`nothing` lets LinearSolve choose).
        """
        struct $name{S<:AbstractTransportSolver,L} <: AbstractAeolisAlgorithm
            solver::S
            linsolve::L
        end
        $name(solver::AbstractTransportSolver=Trunk(); linsolve=nothing) = $name(solver,
                                                                                 linsolve)
        _scheme_symbol(::$name) = $(QuoteNode(sym))
    end
end

"""
    default_algorithm(p::AeolisParameters)

The algorithm named by the `scheme` and `solver` parameters.
"""
function default_algorithm(p::AeolisParameters)
    sol = _solver_from_symbol(p.solver)
    p.scheme == :euler_backward && return EulerBackward(sol)
    p.scheme == :crank_nicolson && return CrankNicolson(sol)
    return EulerForward(sol)
end

"""
    AeolisProblem(p::AeolisParameters; tspan=(p.tstart, p.tstop), u0=nothing)
    AeolisProblem(configfile::AbstractString; kwargs...)

An AeoLiS simulation from `tspan[1]` to `tspan[2]`. Like AeoLiS, steps are taken while the
time is at most `tspan[2]`. `u0` replaces the initial state built from the parameters. A
configuration file is read with [`read_config`](@ref).
"""
struct AeolisProblem{U} <: SciMLBase.AbstractSciMLProblem
    p::AeolisParameters
    tspan::Tuple{Float64,Float64}
    u0::U
end
function AeolisProblem(p::AeolisParameters; tspan=(p.tstart, p.tstop), u0=nothing)
    tspan = (Float64(tspan[1]), Float64(tspan[2]))
    p.tstart == tspan[1] && p.tstop == tspan[2] ||
        (p = remake_parameters(p; tstart=tspan[1], tstop=tspan[2]))
    return AeolisProblem(p, tspan, u0)
end
AeolisProblem(file::AbstractString; kwargs...) = AeolisProblem(read_config(file); kwargs...)

"""
    remake(prob::AeolisProblem; p=prob.p, tspan=prob.tspan, u0=prob.u0, kwargs...)

A copy of `prob` with parts replaced; other keywords replace parameters.
"""
function SciMLBase.remake(prob::AeolisProblem; p=prob.p, tspan=prob.tspan, u0=prob.u0,
                          kwargs...)
    isempty(kwargs) || (p = remake_parameters(p; kwargs...))
    return AeolisProblem(p; tspan, u0)
end

function Base.show(io::IO, ::MIME"text/plain", prob::AeolisProblem)
    println(io, "AeolisProblem with tspan ", prob.tspan)
    return show(io, MIME"text/plain"(), prob.p)
end

# ----------------------------------------------------------------------------------------
# Output

const STATISTICS = (:avg, :sum, :var, :min, :max)

function _parse_output_vars(vars, types)
    out = Dict{Symbol,Vector{Union{Nothing,Symbol}}}()
    order = Symbol[]
    for v in vars
        sv = String(v)
        name, ext = if occursin('_', sv) && Symbol(split(sv, '_')[end]) in STATISTICS
            parts = split(sv, '_')
            Symbol(join(parts[1:(end - 1)], '_')), Symbol(parts[end])
        else
            Symbol(sv), nothing
        end
        haskey(STATE_INDEX, name) || throw(ArgumentError("unknown output variable :$name"))
        haskey(out, name) || (out[name]=Union{Nothing,Symbol}[]; push!(order, name))
        ext in out[name] || push!(out[name], ext)
        for t in types
            t in STATISTICS ||
                throw(ArgumentError("unknown statistic :$t; statistics are $(STATISTICS)"))
            t in out[name] || push!(out[name], t)
        end
    end
    return [(n => out[n]) for n in order]
end

mutable struct _Stats
    min::Array{Float64}
    max::Array{Float64}
    var::Array{Float64}
    sum::Array{Float64}
end
_Stats(sz) = _Stats(fill(Inf, sz), fill(-Inf, sz), zeros(sz), zeros(sz))

"""
    AeolisSolution

The result of solving an [`AeolisProblem`](@ref): output times `t`, and for every saved
variable an array with time as the last dimension. `sol[:zb]` gives the bed level at every
output time, `sol[:zb, k]` at the `k`-th; statistics are named `:zb_avg`, `:zb_max` and so
on. `sol.u` is the final state, `sol.model` the model, `sol.retcode` a SciMLBase `ReturnCode`.
"""
struct AeolisSolution{T,M,P}
    t::Vector{Float64}
    saved::Dict{Symbol,Array{T}}
    u::AeolisState{T}
    model::M
    prob::P
    alg::AbstractAeolisAlgorithm
    retcode::SciMLBase.ReturnCode.T
    stats::Dict{Symbol,Int}
end

Base.getindex(sol::AeolisSolution, name::Symbol) = sol.saved[name]
function Base.getindex(sol::AeolisSolution, name::Symbol, k::Integer)
    a = sol.saved[name]
    return selectdim(a, ndims(a), k)
end
Base.keys(sol::AeolisSolution) = keys(sol.saved)
Base.haskey(sol::AeolisSolution, name::Symbol) = haskey(sol.saved, name)
Base.length(sol::AeolisSolution) = length(sol.t)

function Base.show(io::IO, ::MIME"text/plain", sol::AeolisSolution)
    println(io, "AeolisSolution with ", length(sol.t), " output times from ", first(sol.t),
            " to ", last(sol.t), " s, retcode ", sol.retcode)
    return print(io, "  saved: ", join(sort!(String.(collect(keys(sol.saved)))), ", "))
end
function Base.show(io::IO, sol::AeolisSolution)
    return print(io, "AeolisSolution(", length(sol.t), " times)")
end

"""
    transect(sol, name)

For a transect (a quasi-2-D grid of three rows), the saved variable `name` along the middle
row, as AeoLiS writes it: an array of size `(nx+1, ..., ntimes)`.
"""
function transect(sol::AeolisSolution, name::Symbol)
    a = sol.saved[name]
    size(a, 2) == 3 ||
        throw(ArgumentError("the solution is not a transect (ny + 1 = $(size(a, 2)))"))
    return dropdims(selectdim(a, 2, 2:2); dims=2)
end

# ----------------------------------------------------------------------------------------
# Integrator

"""
    AeolisIntegrator

The state of a running simulation, from `init(prob, alg)`: the state `u` (an
[`AeolisState`](@ref)), the time `t`, the last step `dt`, the model with parameters
`model.p`, and the saved output. `integ[:zb]` reads a state variable; `step!` advances it.
"""
mutable struct AeolisIntegrator{T,M,P,A,C}
    prob::P
    alg::A
    model::M
    u::AeolisState{T}
    t::Float64
    dt::Float64
    dtreq::Float64
    callback::C
    tstops::Vector{Float64}
    saveat::Any
    save_list::Vector{Pair{Symbol,Vector{Union{Nothing,Symbol}}}}
    stats_acc::Dict{Symbol,_Stats}
    nstat::Int
    tout::Float64
    tsaved::Vector{Float64}
    saved::Dict{Symbol,Vector{Array{T}}}
    fixed::Dict{Symbol,Array{T}}
    restart_every::Union{Nothing,Float64}
    restart_file::Union{Nothing,String}
    trestart::Float64
    u_modified::Bool
    terminated::Bool
    retcode::SciMLBase.ReturnCode.T
end

Base.getindex(integ::AeolisIntegrator, name::Symbol) = getfield(integ.u, name)
function Base.setindex!(integ::AeolisIntegrator, v, name::Symbol)
    getfield(integ.u, name) .= v
    integ.u_modified = true
    return integ
end

"""
    statistic(integ, name, stat)

The running statistic `stat` (`:avg`, `:sum`, `:var`, `:min` or `:max`) of state variable
`name` since the last output, as AeoLiS's `get_statistic` (`:sum` is multiplied by the last
time step).
"""
function statistic(integ::AeolisIntegrator, name::Symbol, stat::Symbol)
    o = integ.stats_acc[name]
    n = integ.nstat
    stat == :min && return copy(o.min)
    stat == :max && return copy(o.max)
    stat == :sum && return o.sum .* integ.dt
    stat == :avg && return n > 0 ? o.sum ./ n : zeros(size(o.sum))
    stat == :var &&
        return n > 1 ? (o.var .- o.sum .^ 2 ./ n) ./ (n - 1) : zeros(size(o.var))
    return throw(ArgumentError("unknown statistic :$stat; statistics are $(STATISTICS)"))
end

function _stats_clear!(integ)
    for (name, _) in integ.save_list
        integ.stats_acc[name] = _Stats(size(getfield(integ.u, name)))
    end
    integ.nstat = 0
    return integ
end

function _stats_update!(integ)
    for (name, exts) in integ.save_list
        v = getfield(integ.u, name)
        o = integ.stats_acc[name]
        :min in exts && (o.min .= min.(o.min, v))
        :max in exts && (o.max .= max.(o.max, v))
        any(in(exts), (:sum, :avg, :var)) && (o.sum .+= v)
        :var in exts && (o.var .+= v .^ 2)
    end
    integ.nstat += 1
    return integ
end

function _save!(integ; force::Bool=false)
    t = integ.t
    due = if integ.saveat isa Real
        (t - integ.tout >= integ.saveat) || t == 0
    else
        # a vector of times saves at those times and, as in SciML, at the start
        any(x -> isapprox(x, t; atol=1e-9), integ.saveat) || isempty(integ.tsaved)
    end
    (due || force) || return integ
    push!(integ.tsaved, t)
    for (name, exts) in integ.save_list
        for ext in exts
            key = ext === nothing ? name : Symbol(name, :_, ext)
            v = ext === nothing ? copy(getfield(integ.u, name)) :
                statistic(integ, name, ext)
            push!(get!(integ.saved, key, Array{eltype(integ.u.zb)}[]), v)
        end
    end
    _stats_clear!(integ)
    integ.tout = t
    if integ.restart_every !== nothing && integ.restart_file !== nothing &&
       t - integ.trestart >= integ.restart_every
        save_restart(integ, integ.restart_file)
        integ.trestart = t
    end
    return integ
end

function _discrete_callbacks(cb)
    cb === nothing && return ()
    if cb isa SciMLBase.CallbackSet
        isempty(cb.continuous_callbacks) ||
            throw(ArgumentError("continuous callbacks are not supported; use DiscreteCallback"))
        return cb.discrete_callbacks
    elseif cb isa SciMLBase.DiscreteCallback
        return (cb,)
    end
    return throw(ArgumentError("callback must be a DiscreteCallback or CallbackSet, got $(typeof(cb))"))
end

"""
    init(prob::AeolisProblem, alg=default_algorithm(prob.p); saveat=prob.p.output_times,
         save_vars=prob.p.output_vars, save_stats=prob.p.output_types, callback=nothing,
         tstops=Float64[], dt=-1, fixed=Dict(), restart=nothing, restart_every=prob.p.restart,
         restart_file=nothing, rng=Random.Xoshiro(0), T=Float64)

An [`AeolisIntegrator`](@ref) for `prob`. Output follows AeoLiS: a save whenever at least
`saveat` seconds have passed since the last one (and at `t == 0`); a vector of times saves at
those times, which are added to `tstops`. `save_vars` are state variables, optionally with a
statistic suffix (`:zb_avg`); `save_stats` adds statistics to every saved variable. `dt` is
the requested step (`-1` for the `dt` parameter). Variables in `fixed` keep the given values
(AeoLiS hotstart files). `restart` is a file written by [`save_restart`](@ref).
"""
function CommonSolve.init(prob::AeolisProblem,
                          alg::AbstractAeolisAlgorithm=default_algorithm(prob.p);
                          saveat=prob.p.output_times, save_vars=prob.p.output_vars,
                          save_stats=prob.p.output_types, callback=nothing,
                          tstops=Float64[], dt::Real=-1, fixed=Dict{Symbol,Any}(),
                          restart=nothing, restart_every=prob.p.restart,
                          restart_file=nothing,
                          rng::Random.AbstractRNG=Random.Xoshiro(0),
                          T::Type{<:AbstractFloat}=Float64)
    p = prob.p
    model, s = setup_model(p; T, rng)
    prob.u0 === nothing || copyto_state!(s, prob.u0)
    t0 = prob.tspan[1]
    if restart !== nothing
        t0 = _load_restart!(s, model, restart)
    end
    model.time = t0
    tst = sort!(Float64.(collect(tstops)))
    saveat isa AbstractVector && append!(tst, Float64.(saveat))
    sort!(unique!(tst))
    save_list = _parse_output_vars(Symbol.(collect(save_vars)),
                                   Symbol.(collect(save_stats)))
    fx = Dict{Symbol,Array{T}}()
    for (k, v) in fixed
        a = getfield(s, Symbol(k))
        a .= v
        fx[Symbol(k)] = copy(a)
    end
    integ = AeolisIntegrator(prob, alg, model, s, t0, 0.0, Float64(dt), callback, tst,
                             saveat isa AbstractVector ? Float64.(saveat) : Float64(saveat),
                             save_list, Dict{Symbol,_Stats}(), 0, 0.0, Float64[],
                             Dict{Symbol,Vector{Array{T}}}(), fx,
                             restart_every === nothing ? nothing : Float64(restart_every),
                             restart_file, 0.0, false, false, SciMLBase.ReturnCode.Default)
    _stats_clear!(integ)
    for cb in _discrete_callbacks(callback)
        cb.initialize(cb, integ.u, integ.t, integ)
    end
    _save!(integ)
    return integ
end

SciMLBase.u_modified!(integ::AeolisIntegrator, b::Bool) = (integ.u_modified=b; nothing)
function SciMLBase.derivative_discontinuity!(integ::AeolisIntegrator, b::Bool)
    return SciMLBase.u_modified!(integ, b)
end
function SciMLBase.terminate!(integ::AeolisIntegrator,
                              retcode::SciMLBase.ReturnCode.T=SciMLBase.ReturnCode.Terminated)
    integ.terminated = true
    integ.retcode = retcode
    return nothing
end
function SciMLBase.add_tstop!(integ::AeolisIntegrator, t)
    push!(integ.tstops, Float64(t))
    sort!(unique!(integ.tstops))
    return nothing
end
function SciMLBase.get_proposed_dt(integ::AeolisIntegrator)
    return integ.dtreq > 0 ? integ.dtreq : integ.model.p.dt
end

function _next_step(integ)
    dt = integ.dtreq
    filter!(>(integ.t + 1e-9), integ.tstops)
    if !isempty(integ.tstops)
        step = dt > 0 ? dt : integ.model.p.dt
        integ.t + step > integ.tstops[1] && (dt = integ.tstops[1] - integ.t)
    end
    return dt
end

"""
    step!(integ::AeolisIntegrator, n=1)

Take `n` AeoLiS updates, saving output and applying callbacks after each.
"""
function CommonSolve.step!(integ::AeolisIntegrator, n::Integer=1)
    alg = integ.alg
    for _ in 1:n
        integ.terminated && break
        dt = _next_step(integ)
        step = aeolis_step!(integ.u, integ.model, integ.t, dt, _scheme_symbol(alg),
                            _solver_symbol(alg.solver), alg.linsolve)
        if step === nothing
            integ.t += abs(dt)
            continue
        end
        integ.dt = step
        integ.t += step * integ.model.p.accfac
        for (k, v) in integ.fixed
            getfield(integ.u, k) .= v
        end
        _stats_update!(integ)
        for cb in _discrete_callbacks(integ.callback)
            if cb.condition(integ.u, integ.t, integ)
                cb.affect!(integ)
            end
        end
        _save!(integ)
    end
    return integ
end

"""
    solve!(integ::AeolisIntegrator)

Step while the time is at most the end of the time span, then return the
[`AeolisSolution`](@ref).
"""
function CommonSolve.solve!(integ::AeolisIntegrator)
    tend = integ.prob.tspan[2]
    while integ.t <= tend && !integ.terminated
        step!(integ)
    end
    integ.retcode == SciMLBase.ReturnCode.Default &&
        (integ.retcode = SciMLBase.ReturnCode.Success)
    return _build_solution(integ)
end

function _build_solution(integ::AeolisIntegrator{T}) where {T}
    saved = Dict{Symbol,Array{T}}()
    for (k, v) in integ.saved
        isempty(v) && continue
        saved[k] = cat(v...; dims=ndims(first(v)) + 1)
    end
    return AeolisSolution(copy(integ.tsaved), saved, copy(integ.u), integ.model, integ.prob,
                          integ.alg, integ.retcode, copy(integ.model.counters))
end

"""
    solve(prob::AeolisProblem, alg=default_algorithm(prob.p); kwargs...)

Run the simulation and return an [`AeolisSolution`](@ref); keywords as for `init`.
"""
function CommonSolve.solve(prob::AeolisProblem,
                           alg::AbstractAeolisAlgorithm=default_algorithm(prob.p);
                           kwargs...)
    return solve!(init(prob, alg; kwargs...))
end

"""
    reinit!(integ::AeolisIntegrator)

Restart the integrator from the initial state of its problem, clearing saved output.
"""
function SciMLBase.reinit!(integ::AeolisIntegrator)
    model, s = setup_model(integ.prob.p; T=eltype(integ.u.zb), rng=integ.model.rng)
    integ.prob.u0 === nothing || copyto_state!(s, integ.prob.u0)
    integ.model = model
    integ.u = s
    integ.t = integ.prob.tspan[1]
    integ.dt = 0.0
    integ.tout = 0.0
    integ.terminated = false
    integ.retcode = SciMLBase.ReturnCode.Default
    empty!(integ.tsaved)
    empty!(integ.saved)
    _stats_clear!(integ)
    _save!(integ)
    return integ
end

# ----------------------------------------------------------------------------------------
# Restart files

"""
    save_restart(integ, file)

Write the time, state and counters of `integ` to `file` (Julia serialization), from which
`init(prob; restart=file)` continues.
"""
function save_restart(integ::AeolisIntegrator, file::AbstractString)
    open(file, "w") do io
        return Serialization.serialize(io,
                                       (; t=integ.t, u=integ.u, prev=integ.model.prev,
                                        counters=integ.model.counters,
                                        dt_opt=integ.model.dt_opt))
    end
    return file
end

function _load_restart!(s, model, file)
    data = open(Serialization.deserialize, file)
    copyto_state!(s, data.u)
    copyto_state!(model.prev, data.prev)
    merge!(model.counters, data.counters)
    model.dt_opt = data.dt_opt
    return data.t
end
