function beach_parameters(; kwargs...)
    x = collect(0.0:1.0:60.0)
    zb = @. -1.0 + 0.05 * x
    return AeolisParameters(; xgrid_file=x, bed_file=zb,
                            wind_file=[0.0 12.0 270.0; 86400.0 12.0 270.0],
                            tide_file=[0.0 0.0; 86400.0 0.0], process_tide=true, dt=600.0,
                            tstop=3600.0, output_times=1200.0,
                            output_vars=[:zb, :Ct, :pickup],
                            kwargs...)
end

@testset "problem interface" begin
    p = beach_parameters()
    prob = AeolisProblem(p)
    @test prob.tspan == (0.0, 3600.0)
    @test default_algorithm(p) isa EulerBackward{Trunk}
    sol = solve(prob)
    @test sol.retcode == ReturnCode.Success
    # AeoLiS saves at t = 0 and whenever output_times has passed; it steps while t <= tstop
    @test sol.t == [0.0, 1200.0, 2400.0, 3600.0]
    @test size(sol[:zb]) == (61, 3, 4)
    @test size(sol[:Ct]) == (61, 3, 1, 4)
    @test sol[:zb, 1] == p.bed_file .* ones(1, 3)
    @test size(transect(sol, :zb)) == (61, 4)
    @test sum(sol[:pickup]) > 0
    @test occursin("4 output times", sprint(show, MIME"text/plain"(), sol))
    # the same with an explicit algorithm and LinearSolve's choice
    sol2 = solve(prob, EulerBackward(Trunk()))
    @test sol2[:zb] == sol[:zb]
    # statistics
    sol3 = solve(prob; save_vars=[:zb, :Ct], save_stats=[:avg, :max])
    @test haskey(sol3, :Ct_avg) && haskey(sol3, :Ct_max)
    # at the first output no step has been taken: AeoLiS's max starts at -Inf, the average at 0
    @test all(sol3[:Ct_max][:, :, :, 2:end] .>= sol3[:Ct_avg][:, :, :, 2:end] .- 1e-15)
    # init / step!
    integ = init(prob)
    @test integ.t == 0.0
    step!(integ)
    @test integ.t == 600.0
    @test integ.dt == 600.0
    step!(integ, 2)
    @test integ.t == 1800.0
    @test integ[:zb] === integ.u.zb
    zb = copy(integ[:zb])
    integ[:zb] = zb .+ 0.01
    @test integ[:zb] ≈ zb .+ 0.01
    # remake
    prob2 = remake(prob; tspan=(0.0, 1200.0), grain_size=0.0004)
    @test prob2.p.grain_size == [0.0004]
    @test solve(prob2).t == [0.0, 1200.0]
    # callbacks: a DiscreteCallback that terminates, and one that changes the state
    hits = Float64[]
    cb = DiscreteCallback((u, t, integ) -> t >= 1800,
                          integ -> (push!(hits, integ.t); terminate!(integ)))
    solc = solve(prob; callback=cb)
    @test solc.retcode == ReturnCode.Terminated
    @test hits == [1800.0]
    raise = DiscreteCallback((u, t, integ) -> t == 600,
                             integ -> (integ[:zb] = integ[:zb] .+ 1.0))
    solr = solve(prob; callback=CallbackSet(raise))
    @test maximum(abs.(solr[:zb, 2] .- (sol[:zb, 2] .+ 1.0))) < 0.01
    @test_throws ArgumentError solve(prob; callback=1)
    # tstops shorten the step to land on them
    integ = init(prob; tstops=[900.0])
    step!(integ, 2)
    @test integ.t == 900.0
    # saveat as a vector of times
    solv = solve(prob; saveat=[1500.0, 3000.0])
    @test solv.t == [0.0, 1500.0, 3000.0]
    # fixed variables (AeoLiS hotstart files) are kept
    solf = solve(prob; fixed=Dict(:zb => p.bed_file .* ones(1, 3)))
    @test all(solf[:zb, k] == solf[:zb, 1] for k in 2:4)
    # restart files
    integ = init(prob)
    step!(integ, 3)
    file = joinpath(mktempdir(), "restart.jls")
    save_restart(integ, file)
    rest = solve(prob; restart=file)
    @test rest.t[1] == 1800.0
    @test rest.u.zb ≈ sol.u.zb
    # reinit!
    integ = init(prob)
    step!(integ, 2)
    reinit!(integ)
    @test integ.t == 0.0 && integ.u.zb == sol[:zb, 1]
    # ensembles through SciMLBase
    eprob = EnsembleProblem(prob;
                            prob_func=(pr, ctx) -> remake(pr;
                                                          grain_size=[0.0002, 0.0004][ctx.sim_id]))
    esol = solve(eprob, EulerBackward(), Aeolis.SciMLBase.EnsembleSerial(); trajectories=2)
    @test length(esol.u) == 2
    @test sum(esol.u[1][:pickup]) > sum(esol.u[2][:pickup])
    # Float32 states
    sol32 = solve(prob; T=Float32)
    @test eltype(sol32[:zb]) == Float32
    @test sol32[:zb] ≈ sol[:zb] rtol = 1e-4
end

@testset "schemes and solvers" begin
    p = beach_parameters()
    prob = AeolisProblem(p)
    ref = solve(prob)[:zb]
    for alg in (CrankNicolson(), EulerBackward(Pieter()), EulerBackward(SteadyState()),
                EulerForward(Trunk()))
        s = solve(remake(prob; CFL=1.0), alg)
        @test s.retcode == ReturnCode.Success
        @test maximum(abs.(s[:zb][:, :, end] .- ref[:, :, end])) < 0.05
    end
    @test_throws UnsupportedOptionError solve(prob, CrankNicolson(SteadyState()))
end
