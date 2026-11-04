@testset "netCDF output" begin
    x = collect(0.0:1.0:30.0)
    p = AeolisParameters(; xgrid_file=x, bed_file=0.05 .* x,
                         wind_file=[0.0 12.0 270.0; 7200.0 12.0 270.0], dt=600.0,
                         tstop=1200.0, output_times=600.0,
                         output_vars=[:zb, :Ct, :mass], output_types=[:avg])
    sol = solve(AeolisProblem(p))
    file = joinpath(mktempdir(), "out.nc")
    write_netcdf(file, sol)
    NCDataset(file) do ds
        @test Array(ds["time"].var) == sol.t
        @test size(ds["zb"]) == (31, 1, 4)          # transects as one row; steps while t <= tstop
        @test ds["zb"][:, 1, :] ≈ Float32.(transect(sol, :zb))
        @test size(ds["mass"]) == (1, 3, 31, 1, 4)   # (fractions, layers, s, n, time)
        @test haskey(ds, "zb_avg")
        @test ds["zb"].attrib["units"] == "m"
    end
end
