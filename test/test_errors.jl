@testset "errors" begin
    @test_throws ConfigurationError solve(AeolisProblem(AeolisParameters()))
    x = collect(0.0:1.0:10.0)
    base = (; xgrid_file=x, bed_file=zeros(11))
    @test_throws ConfigurationError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        wind_file=[0.0 1.0]))
    @test_throws UnknownOptionError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        scheme=:leapfrog))
    @test_throws UnknownOptionError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        method_transport=:magic))
    @test_throws UnsupportedOptionError Aeolis.setup_model(AeolisParameters(; base...,
                                                                            th_drylayer=true))
    @test_throws UnsupportedOptionError Aeolis.setup_model(AeolisParameters(; base...,
                                                                            th_humidity=true))
    @test_throws ConfigurationError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        bed_file=zeros(5)))
    @test_throws ConfigurationError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        solver=:steadystate,
                                                                        grain_size=[1e-4,
                                                                                    2e-4],
                                                                        grain_dist=[0.5,
                                                                                    0.5]))
    e = UnknownOptionError(:scheme, :leapfrog, [:euler_forward])
    @test occursin(":leapfrog", sprint(showerror, e)) &&
          occursin(":euler_forward", sprint(showerror, e))
    e = UnsupportedOptionError(:solver, :x, "reason")
    @test occursin("reason", sprint(showerror, e))
    @test occursin("ConfigurationError: msg", sprint(showerror, ConfigurationError("msg")))
    @test ConfigurationError <: AeolisError
    @test_throws ArgumentError solve(AeolisProblem(AeolisParameters(; base...,
                                                                    wind_file=[0.0 5.0 270.0;
                                                                               3600.0 5.0 270.0]));
                                     save_vars=[:nonsense])
end
