@testset "parameters" begin
    p = AeolisParameters()
    @test p.grain_size == [225e-6]
    @test p.scheme == :euler_backward
    @test p.output_vars == [:zb, :zs, :Ct, :Cu, :uw, :udir, :uth, :mass, :pickup, :w]
    @test p.Cs == 2.5e-3
    q = AeolisParameters(; scheme="Crank_Nicolson", grain_size=0.0003, nlayers=2.0,
                         output_vars="zb", process_tide="T")
    @test q.scheme == :crank_nicolson
    @test q.grain_size == [0.0003]
    @test q.nlayers == 2
    @test q.output_vars == [:zb]
    @test q.process_tide
    @test_throws ArgumentError AeolisParameters(; grain_szie=1.0)
    err = try
        AeolisParameters(; grain_szie=1.0)
    catch e
        e
    end
    @test occursin("grain_size", sprint(showerror, err))
    @test_throws ArgumentError AeolisParameters(; nlayers="many")
    r = remake_parameters(q; tstop=10.0)
    @test r.tstop == 10.0 && r.scheme == :crank_nicolson
    @test q == remake_parameters(q)
    @test length(parameter_info()) == length(fieldnames(AeolisParameters))
    @test parameter_info(:Cb).default == 1.5
    @test_throws ArgumentError parameter_info(:nope)
    @test state_info(:zb).unit == "m"
    @test_throws ArgumentError state_info(:nope)
    s = AeolisState{Float64}(4, 2, 3, 2)
    @test size(s.zb) == (5, 3)
    @test size(s.Ct) == (5, 3, 2)
    @test size(s.mass) == (5, 3, 3, 2)
    s[:zb] = 1.0
    @test all(==(1.0), s[:zb])
    @test sprint(show, MIME"text/plain"(), s) ==
          "AeolisState{Float64} on 5 x 3 points, 3 layers, 2 fractions"
    @test occursin("scheme: crank_nicolson", sprint(show, MIME"text/plain"(), q))
end

@testset "configuration files" begin
    dir = mktempdir()
    x = collect(0.0:1.0:10.0)
    open(joinpath(dir, "x.grd"), "w") do io
        println(io, join(x, " "))
    end
    open(joinpath(dir, "z.grd"), "w") do io
        println(io, join(0.1 .* x, " "))
    end
    open(joinpath(dir, "wind.txt"), "w") do io
        println(io, "0 10 270")
        println(io, "3600 10 270")
    end
    cfg = joinpath(dir, "aeolis.txt")
    open(cfg, "w") do io
        println(io, "% a comment line = ignored")
        println(io, "xgrid_file = x.grd    % grid")
        println(io, "bed_file = z.grd")
        println(io, "wind_file = wind.txt")
        println(io, "grain_size = 0.0002 0.0004")
        println(io, "grain_dist = 0.5 0.5")
        println(io, "max_error = 1e-7")
        println(io, "refdate = 2020-01-01 00:00")
        println(io, "process_tide = T")
        println(io, "Tswash = 3")
    end
    p = @test_logs (:warn,) read_config(cfg)
    @test p.xgrid_file == x
    @test size(p.wind_file) == (2, 3)
    @test p.grain_size == [0.0002, 0.0004]
    @test p.max_error == 1e-7
    @test p.refdate == "2020-01-01 00:00"
    @test p.process_tide
    @test p.output_file == joinpath(dir, "aeolis.nc")
    @test p.h == p.z
    p2 = read_config(cfg; tstop=60.0)
    @test p2.tstop == 60.0
    out = joinpath(dir, "copy", "aeolis.txt")
    mkpath(dirname(out))
    write_config(out, p)
    p3 = read_config(out)
    @test p3.xgrid_file ≈ p.xgrid_file
    @test p3.grain_size == p.grain_size
    @test p3.wind_file ≈ p.wind_file
    @test_throws ConfigurationError read_config(joinpath(dir, "missing.txt"))
end
