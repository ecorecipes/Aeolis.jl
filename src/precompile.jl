# A small workload run when the package is precompiled, so that the first model run does not
# spend most of a minute compiling the time step: a transect with wind, tide and avalanching,
# a transect with marram vegetation, and a small 2-D grid with spectral shear.
PrecompileTools.@setup_workload begin
    x = collect(0.0:1.0:40.0)
    zb = @. -1.0 + 0.1 * x
    wind = [0.0 12.0 270.0; 86400.0 12.0 270.0]
    tide = [0.0 0.0; 86400.0 0.0]
    X = [i for i in 0.0:2.0:30.0, j in 0.0:2.0:20.0]
    Y = [j for i in 0.0:2.0:30.0, j in 0.0:2.0:20.0]
    Z = @. 2.0 * exp(-((X - 15)^2 + (Y - 10)^2) / 30)
    PrecompileTools.@compile_workload begin
        solve(AeolisProblem(AeolisParameters(; xgrid_file=x, bed_file=zb, wind_file=wind,
                                             tide_file=tide, process_tide=true,
                                             process_avalanche=true, dt=600.0, tstop=1800.0,
                                             output_times=600.0)))
        solve(AeolisProblem(AeolisParameters(; xgrid_file=x, bed_file=zb, wind_file=wind,
                                             tide_file=tide, process_tide=true,
                                             process_vegetation=true, veg_model=:marram,
                                             veg_file=(@. ifelse(30 <= x <= 35, 0.3, 0.0)),
                                             dt=3600.0, tstop=7200.0, output_times=3600.0)))
        solve(AeolisProblem(AeolisParameters(; xgrid_file=X, ygrid_file=Y, bed_file=Z,
                                             wind_file=wind, process_shear=true,
                                             process_avalanche=true, dt=600.0, tstop=1200.0,
                                             output_times=600.0)))
    end
end
