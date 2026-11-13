# SciMLSteadyState: the equilibrium concentration as a SteadyStateProblem.

function flat_parameters(; kwargs...)
    x = collect(0.0:1.0:40.0)
    return AeolisParameters(; xgrid_file=x, bed_file=fill(5.0, length(x)),
                            wind_file=[0.0 10.0 0.0; 1e6 10.0 0.0], tstart=0.0,
                            tstop=1200.0, dt=600.0, output_times=600.0,
                            wind_convention=:cartesian, boundary_offshore=:flux,
                            offshore_flux=0.2, boundary_onshore=:gradient,
                            boundary_lateral=:circular, kwargs...)
end

# net pickup per second against outflow minus inflow through the cross-shore boundaries
function mass_balance(s, p, dt)
    C = s.Ct[:, :, 1]
    us = s.us[:, :, 1]
    pick = sum(s.pickup[:, :, 1] .* s.dsdn) / dt
    Cin = p.offshore_flux .* s.Cu0[1, :, 1]
    inflow = sum(s.dn[1, :] .* max.(us[1, :], 0) .* Cin) +
             sum(s.dn[1, :] .* max.(-us[1, :], 0) .* C[1, :])
    outflow = sum(s.dn[end, :] .* us[end, :] .* C[end, :])
    return pick, outflow - inflow
end

@testset "SciMLSteadyState" begin
    alg = EulerBackward(SciMLSteadyState(NewtonRaphson(); abstol=1e-14))
    @test alg.solver.kwargs == (; abstol=1e-14)

    @testset "exact upwind solution on a flat transect" begin
        p = flat_parameters()
        integ = init(AeolisProblem(p), alg)
        step!(integ, 2)
        s = integ.u
        C = s.Ct[:, 2, 1]
        v = s.us[1, 2, 1]
        @test all(≈(v), s.us[:, 2, 1])
        prev = p.offshore_flux * s.Cu0[1, 2, 1]
        exact = similar(C)
        for a in eachindex(C)
            r = s.ds[a, 2] / (v * p.T)
            exact[a] = (prev + r * s.Cu[a, 2, 1]) / (1 + r)
            prev = exact[a]
        end
        @test maximum(abs.(C .- exact)) <= 1e-10 * maximum(exact)
        pick, net = mass_balance(s, p, integ.model.dt_opt)
        @test pick ≈ net rtol = 1e-8
    end

    @testset "supply limitation" begin
        p = flat_parameters(; layer_thickness=1e-6, nlayers=3)
        integ = init(AeolisProblem(p), alg)
        mass0 = copy(integ.u.mass[:, :, 1, 1])
        step!(integ, 2)
        s = integ.u
        @test integ.model.counters[:supplylim] > 0
        @test all(s.pickup[:, :, 1] .<= mass0 .+ 1e-12)
        @test any(s.w[:, :, 1] .< 1 - 1e-6)
        pick, net = mass_balance(s, p, integ.model.dt_opt)
        @test pick ≈ net rtol = 1e-8
    end

    @testset "two dimensions and fractions" begin
        x, y = 0.0:1.0:30.0, 0.0:1.0:20.0
        X = [xi for xi in x, _ in y]
        Y = [yi for _ in x, yi in y]
        Z = @. 2.0 * exp(-((X - 12)^2 / 30 + (Y - 10)^2 / 40))
        p = AeolisParameters(; xgrid_file=X, ygrid_file=Y, bed_file=Z,
                             wind_file=[0.0 13.0 260.0; 1e6 13.0 260.0], tstop=1200.0,
                             dt=600.0, output_times=600.0, process_shear=true,
                             boundary_offshore=:flux, offshore_flux=0.3,
                             grain_size=[0.0002, 0.0004], grain_dist=[0.5, 0.5])
        integ = init(AeolisProblem(p), alg)
        step!(integ, 2)
        @test all(isfinite, integ.u.Ct)
        @test all(>=(0), integ.u.Ct)
        @test maximum(integ.u.Ct) > 0
    end

    @testset "unsupported boundaries" begin
        p = flat_parameters(; boundary_offshore=:uniform)
        @test_throws UnsupportedOptionError solve(AeolisProblem(p), alg)
    end
end
