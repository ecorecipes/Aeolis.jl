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
