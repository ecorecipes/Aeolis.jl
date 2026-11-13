@testset "process formulas" begin
    # Stockdon et al. (2006): dissipative branch (iribarren < 0.3) and the general one
    eta, sig, R = Aeolis.stockdon_runup(1.0, 8.0, 0.01)
    Lo = 9.81 * 64 / (2pi)
    @test R ≈ 0.043 * sqrt(Lo)
    @test eta ≈ R - sig
    eta, sig, R = Aeolis.stockdon_runup(1.0, 8.0, 0.1)
    @test eta ≈ 0.35 * 0.1 * sqrt(Lo)
    @test R ≈ 1.1 * (eta + sig)
    @test Aeolis.stockdon_runup(0.0, 8.0, 0.1) == (0.0, 0.0, 0.0)
    # a flat transect under a steady wind
    x = collect(0.0:1.0:50.0)
    p = AeolisParameters(; xgrid_file=x, bed_file=fill(1.0, length(x)),
                         wind_file=[0.0 10.0 270.0; 7200.0 10.0 270.0], dt=600.0,
                         tstop=1800.0, method_transport=:bagnold,
                         method_grainspeed=:windspeed)
    model, s = Aeolis.setup_model(p)
    @test size(s.zb) == (51, 3)
    @test s.ds[1] == 1.0 && s.dn[1] == 1.0
    Aeolis.interpolate_wind!(s, model, 0.0)
    @test all(≈(10.0), s.uw)
    @test all(>(0), s.uws)            # a nautical 270 deg wind blows towards +x
    @test s.ustar[1] ≈ 10.0 * 0.41 / log(10.0 / 0.001)
    @test s.tau[1] ≈ p.rhoa * s.ustar[1]^2
    Aeolis.stress_velocity!(s, model)
    @test s.ustar[1] ≈ 10.0 * 0.41 / log(10.0 / 0.001)
    Aeolis.compute_threshold!(s, model)
    @test s.uth[1] ≈ p.Aa * sqrt((p.rhog - p.rhoa) / p.rhoa * p.g * p.grain_size[1])
    Aeolis.equilibrium!(s, model)
    @test s.Cu[1] ≈ p.Cb * p.rhoa / p.g * (s.ustar[1] - s.uth[1])^3 / s.uw[1]
    # bed mass from the grain size distribution
    @test s.mass[1, 1, 1, 1] ≈ p.rhog * (1 - p.porosity) * p.layer_thickness
end

@testset "spectral shear with a spatially varying z0" begin
    @test Aeolis._effective_z0(0.001) === 0.001
    @test Aeolis._effective_z0(fill(2e-5, 3, 4)) === 2e-5
    @test Aeolis._effective_z0([1e-5 NaN; 1e-3 0.0]) ≈ 1e-4
    @test_throws ConfigurationError Aeolis._effective_z0([0.0 NaN])
    # a hill with two fractions: the adaptive roughness is uniform at first, so the first
    # step matches the initial mean grain size (AeoLiS raises a ValueError for the grid);
    # in the second step (AeoLiS steps while t <= tstop) the bed has sorted and the
    # undisturbed shear velocity uses the new z0
    x, y = 0.0:1.0:40.0, 0.0:1.0:30.0
    X = [xi for xi in x, _ in y]
    Y = [yi for _ in x, yi in y]
    Z = @. 3.0 * exp(-((X - 15)^2 / 50 + (Y - 15)^2 / 70))
    hill(meth) = AeolisParameters(; xgrid_file=X, ygrid_file=Y, bed_file=Z,
                                  wind_file=[0.0 13.0 255.0; 1e6 13.0 255.0], tstop=600.0,
                                  dt=600.0, output_times=600.0, process_shear=true, L=50.0,
                                  l=1.0, buffer_width=10.0, grain_size=[0.0002, 0.0004],
                                  grain_dist=[0.5, 0.5], method_roughness=meth)
    adaptive = solve(AeolisProblem(hill(:mean_grainsize_adaptive)))
    initial = solve(AeolisProblem(hill(:mean_grainsize_initial)))
    @test adaptive.retcode == ReturnCode.Success
    @test adaptive[:zb][:, :, 1:2] ≈ initial[:zb][:, :, 1:2] rtol = 1e-10
    @test !(adaptive[:zb][:, :, 3] ≈ initial[:zb][:, :, 3])
    integ = init(AeolisProblem(hill(:mean_grainsize_adaptive)))
    @test integ.model.shear.z0 ≈ (0.0002 + 0.0004) / 2 / 30
end
