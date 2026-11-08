using Dates

@testset "marram traits" begin
    tr = MarramTraits()
    @test tr.N_max == 470.0
    @test length(marram_info()) == fieldcount(MarramTraits)
    @test all(r -> r.evidence in (:M, :D, :P), marram_info())
    @test all(r -> r.range[1] <= r.default <= r.range[2], marram_info())
    @test marram_info(:burial_opt).source == "Nolet2018"
    @test_throws ArgumentError marram_info(:burial_optimum)
    err = try
        MarramTraits(; burial_optimum=0.3)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("burial_opt", sprint(showerror, err))
    @test_throws ArgumentError MarramTraits(; s_w=0.2)
    @test_throws ArgumentError MarramTraits(; tau_P=0.0)
    @test_throws ArgumentError MarramTraits(; r="fast")
    t2 = MarramTraits(tr; r=4.0)
    @test t2.r == 4.0 && t2.N_max == tr.N_max && t2 != tr && MarramTraits(tr) == tr
    @test hash(MarramTraits(tr)) == hash(tr)
    @test occursin("all defaults", sprint(show, MIME"text/plain"(), tr))
    @test occursin("r = 4.0", sprint(show, MIME"text/plain"(), t2))
    mktempdir() do dir
        f = joinpath(dir, "traits.toml")
        write(f, "N_max = 400\n[burial]\nburial_opt = 0.35\n")
        t3 = marram_traits(f)
        @test t3.N_max == 400 && t3.burial_opt == 0.35 && t3.r == tr.r
        write(joinpath(dir, "aeolis.txt"),
              "veg_model = marram\nveg_traits_file = traits.toml\n")
        p = read_config(joinpath(dir, "aeolis.txt"))
        @test p.veg_model == :marram
        @test p.veg_traits_file == normpath(joinpath(dir, "traits.toml"))
        @test marram_traits(p) == t3
        write(f, "N_max = 400\n[shoots]\nN_max = 300\n")
        @test_throws ConfigurationError marram_traits(f)
    end
    @test_throws ConfigurationError marram_traits("no_such_file.toml")
    @test marram_traits(AeolisParameters()) == tr
end

@testset "marram parameters and calendar" begin
    p = AeolisParameters()
    @test p.veg_model == :legacy && p.veg_interval == 86400.0
    x = collect(0.0:1.0:10.0)
    base = (; xgrid_file=x, bed_file=zeros(11),
            wind_file=[0.0 5.0 270.0; 3600.0 5.0 270.0])
    @test_throws UnknownOptionError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        veg_model=:dubeveg))
    @test_throws ConfigurationError Aeolis.setup_model(AeolisParameters(; base...,
                                                                        veg_interval=0.0))
    q = AeolisParameters(; refdate="2016-06-24 12:00")
    @test Aeolis.model_datetime(q, 86400.0) == DateTime(2016, 6, 25, 12)
    @test Aeolis.model_datetime(AeolisParameters(; refdate="2020-03-01"), 0) ==
          DateTime(2020, 3, 1)
    @test_throws ConfigurationError Aeolis.model_datetime(AeolisParameters(;
                                                                           refdate="soon"),
                                                          0)
    @test Aeolis.day_of_year(DateTime(2021, 1, 1, 12)) == 1.5
    @test Aeolis.day_of_year(Date(2021, 4, 1)) == 91.0
end

@testset "marram process functions" begin
    tr = MarramTraits()
    A = Aeolis
    # season: zero outside, ramps, one in the middle
    @test A.season(90.9, tr) == 0 && A.season(274.0, tr) == 0
    @test A.season(91 + 7, tr) ≈ 0.5
    @test A.season(180.0, tr) == 1
    @test A.season(274 - 7, tr) ≈ 0.5
    @test A.season_length(tr) == 183
    @test A.seasonal_burial(0.01, tr) ≈ 1.83
    # burial response: peak at the optimum, split widths
    @test A.burial_response(tr.burial_opt, tr) == 1
    @test A.burial_response(tr.burial_opt - 0.2, tr) ≈ exp(-0.04 / (2 * tr.burial_c_lo^2))
    @test A.burial_response(tr.burial_opt + 0.2, tr) ≈ exp(-0.04 / (2 * tr.burial_c_hi^2))
    @test A.elongation(-0.5, tr) == tr.e0
    # water: continuous, bounded
    for s in (tr.s_w, tr.s_star)
        @test A.water_stress(s - 1e-9, tr) ≈ A.water_stress(s + 1e-9, tr) atol = 1e-6
    end
    @test A.water_stress(0.0, tr) == 0 && A.water_stress(1.0, tr) == 1
    @test all(0 <= A.water_stress(s, tr) <= 1 for s in 0:0.001:1)
    for s in (tr.s_h, tr.s_w, tr.s_star)
        @test A.evapotranspiration(s - 1e-9, 4.0, 0.5, tr) ≈
              A.evapotranspiration(s + 1e-9, 4.0, 0.5, tr) atol = 1e-6
    end
    @test A.evapotranspiration(0.5, 4.0, 1.0, tr) == 4.0
    @test A.evapotranspiration(0.5, 9.0, 1.0, tr) == tr.E_max
    @test A.evapotranspiration(0.5, 4.0, 0.0, tr) == tr.E_w
    @test A.drainage(tr.s_fc, tr) == 0 && A.drainage(1.0, tr) ≈ tr.K_s
    for (rain, dt) in ((0.0, 1.0), (500.0, 1.0), (0.0, 1e4), (1e5, 0.1))
        s = A.bucket_step(0.5, rain, 5.0, 1.0, dt, tr)
        @test tr.s_h <= s <= 1
    end
    @test A.bucket_step(tr.s_fc, 0.0, 0.0, 0.0, 1e4, tr) < tr.s_w   # bare sand dries slowly
    # flooding: survival 1/2 at D50; the hazard integrates to -log survival
    @test A.flood_survival(tr.D50, tr) ≈ 0.5
    @test A.flood_survival(0.0, tr) == 1 && A.flood_hazard(0.0, tr) == 0
    D, h = 10.0, 1e-3
    @test sum(A.flood_hazard((k - 0.5) * h, tr) * h for k in 1:round(Int, D / h)) ≈
          -log(A.flood_survival(D, tr)) rtol = 1e-4
    @test A.drought_hazard(tr.t_dry, tr) == 0
    @test A.drought_hazard(tr.t_dry + 1, tr) == tr.k_dry
    # pressure: the exact solution of its equation, bounded, relieved by burial
    P, d, dt, n = 0.1, 0.0, 400.0, 400_000
    Pe = P
    for _ in 1:n
        Pe += (Pe + tr.P_seed) * (A.pressure_target(d, tr) - Pe) / (tr.tau_P * 365.25) *
              (dt / n)
    end
    @test A.pressure_step(P, d, dt, tr) ≈ Pe rtol = 1e-4
    @test A.pressure_step(0.9, 1.0, 1e5, tr) ≈ 0 atol = 1e-12
    @test A.pressure_step(0.9, 1.0, 365.25 * tr.tau_relief, tr) ≈ 0.9 / ℯ
    @test A.pressure_step(0.0, 0.0, 1e6, tr) ≈ 1
    @test A.pressure_step(0.0, 0.0, 1e300, tr) ≈ 1
    @test A.vigour(1.0, tr) ≈ 1 - tr.p_max
    # roots, reinforcement, drag
    @test A.root_shoot_ratio(0.0, tr) == tr.rs0
    @test A.root_shoot_ratio(tr.rs_d1, tr) ≈ tr.rs1
    @test A.root_shoot_ratio(10.0, tr) == tr.rs2
    @test A.root_target(1e6, 0.5, tr) == tr.Br_max
    @test A.reinforcement(0.0, 0.5, tr) == 0 && A.reinforcement(1e3, 0.5, tr) == tr.c_max
    @test A.reinforcement(0.01, 0.0, tr) == A.reinforcement(0.01, tr.zr_min_cr, tr)
    @test A.standing_fraction(100.0, 0.0) == 0 && A.standing_fraction(0.0, 1e4) == 1
    @test A.frontal_area(0.0, 1.0, 5.0, tr) == 0
    @test A.effective_cover(1.0, tr, 16.0) == 1
    @test A.effective_cover(0.1, tr, 16.0) ≈ tr.m_beta * 0.1 / 16
end

@testset "marram step" begin
    tr = MarramTraits()
    A = Aeolis
    f = (; doy=180.0, dzb=0.0, s=tr.s_fc)
    # without losses, shoot density follows the logistic exactly
    t0 = MarramTraits(tr; m0=0.0, p_max=0.0)
    c = MarramCell(; N=10.0, h=tr.h_max, s=tr.s_fc)
    g = A.marram_rates(c, f, t0).growth
    for dt in (1.0, 100.0)
        N = marram_step(c, f, t0, dt).N
        @test N ≈ tr.N_max / (1 + (tr.N_max / 10 - 1) * exp(-g * dt)) rtol = 1e-10
    end
    # positivity and bounds for any step
    for dt in (1e-3, 1.0, 1e4), dzb in (-0.05, 0.0, 0.05)
        c2 = marram_step(MarramCell(; N=100.0, h=0.5, Br=0.2, zr=0.5),
                         (; doy=180.0, dzb, rain=0.0, pet=5.0, submerged=true), tr, dt)
        @test all(x -> isfinite(x) && x >= 0, (c2.N, c2.h, c2.Br, c2.zr, c2.s, c2.P))
        @test c2.N <= tr.N_max && c2.h <= tr.h_max && c2.P <= 1 && c2.s <= 1
        @test c2.zr <= tr.zr_max && c2.Br <= tr.Br_max
    end
    # erosion of the whole rooted layer removes the plants
    c3 = marram_step(MarramCell(; N=100.0, h=0.5, Br=0.2, zr=0.1),
                     (; doy=180.0, dzb=-0.2), tr, 1.0)
    @test c3.N == 0 && c3.h == 0 && c3.Br == 0
    # burial of the canopy smothers shoots
    r = A.marram_rates(MarramCell(; N=100.0, h=0.0), (; doy=10.0, dzb=0.0), tr)
    @test r.m_sm == 1 / tr.tau_sm && r.growth == 0
    @test marram_rates(c, f, tr).mortality ≈ tr.m0 / 365.25
    # without rain or prescribed moisture, moisture keeps its value
    @test marram_step(MarramCell(; N=1.0, s=0.3), (; doy=1.0, dzb=0.0), tr, 1.0).s == 0.3
    sol = simulate_marram(MarramCell(; N=5.0, h=0.5), (t, d) -> (; dzb=1e-3), tr; days=10)
    @test length(sol.t) == 11 && sol.t[end] == 10
    @test sol.date[end] == DateTime(2020, 1, 11)
    @test sol.N[1] == 5.0 && all(==(0.0), sol.phi)
    @test_throws ArgumentError simulate_marram(MarramCell(), (t, d) -> (; dzb=0.0), tr;
                                               days=0)
end

@testset "marram experiments" begin
    @test_throws ArgumentError marram_experiment(:nope)
    e = marram_experiment(:reijers2021)
    @test all(abs.(e.shoots ./ e.shoots_observed .- 1) .< 0.1)
    @test all(abs.(e.ratio ./ e.ratio_observed .- 1) .< 0.2)
    @test all(abs.(e.elongation .- e.elongation_observed) .< 0.001)
    e = marram_experiment(:nolet2018)
    below = e.burial .<= MarramTraits().burial_opt
    ndvi = vec(sum(e.ndvi_observed; dims=2)) ./ 4
    @test sqrt(sum(abs2, e.response[below] .- ndvi[below]) / count(below)) < 0.03
    @test maximum(abs.(e.response[below] .- ndvi[below])) < 0.08
    lo, hi = vec(minimum(e.ndvi_observed; dims=2)), vec(maximum(e.ndvi_observed; dims=2))
    @test all(lo .- 0.05 .<= e.response .<= hi .+ 0.05)
    e = marram_experiment(:dixon2004)
    @test abs(e.survival - e.survival_observed) < 0.05
    e = marram_experiment(:aptekar2000)
    @test all(abs.(e.survival[3:4] .- e.survival_observed[3:4]) .< 0.02)
    e = marram_experiment(:lammers2024)
    @test e.growth_ratio > 0.95 && e.survival == 1
    e = marram_experiment(:homberger2024)
    @test e.shoots[1] < 0.5 && all(e.shoots[2:end] .> 1)
end

@testset "marram grid model" begin
    day = 86400.0
    x = collect(0.0:0.5:30.0)
    calm = [0.0 2.0 270.0; 4000day 2.0 270.0]       # below the transport threshold
    base = (; xgrid_file=x, bed_file=fill(3.0, length(x)), wind_file=calm, dt=day,
            process_vegetation=true, veg_model=:marram, veg_file=Float64.(x .<= 10),
            refdate="2020-01-01 00:00")
    tr = MarramTraits()
    m, s = Aeolis.setup_model(AeolisParameters(; base..., tstop=10day))
    @test m.veg.traits == tr
    @test s.veg_N[1, 2] == tr.N_max && s.veg_N[end, 2] == 0
    @test s.hveg[1, 2] == tr.h_new && s.veg_s[1, 2] == tr.s_fc
    @test all(>(0), s.veg_E) && all(s.veg_E[:, 1] .== s.veg_E[:, 2] .== s.veg_E[:, 3])
    @test 0 < s.rhoveg[1, 2] <= 1 && s.rhoveg[end, 2] == 0
    # the patch spreads at about 0.4 m per year (Keijsers et al. 2015)
    p = AeolisParameters(; base..., tstop=8 * 365day, output_times=365day,
                         output_vars=[:veg_N], veg_random=false)
    sol = solve(AeolisProblem(p))
    N = transect(sol, :veg_N)
    front = [maximum(x[N[:, k] .> 0]) for k in axes(N, 2)]
    @test 0.3 <= (front[end] - front[3]) / 6 <= 0.5
    integ = init(AeolisProblem(AeolisParameters(; base..., tstop=200day)))
    solve!(integ)
    u = integ.u
    for v in (:veg_N, :hveg, :veg_Br, :veg_P, :rhoveg)
        a = getfield(u, v)
        @test a[:, 1] == a[:, 2] == a[:, 3]
    end
    @test all(isfinite, u.veg_N) && all(0 .<= u.veg_P .<= 1)
    # seawater kills plants that stay under water, and only those
    tide = [0.0 2.5; 4000day 2.5]
    zb = collect(range(2.0, 3.0; length=length(x)))
    q = AeolisParameters(; base..., bed_file=zb, veg_file=ones(length(x)),
                         process_tide=true,
                         tide_file=tide, tstop=40day)
    integ = init(AeolisProblem(q))
    solve!(integ)
    @test all(==(0), integ.u.veg_N[zb .< 2.45, 2])
    @test all(>(0), integ.u.veg_N[zb .> 2.55, 2])
    # without rain, a closed canopy dries its root zone and stops growing
    w = AeolisParameters(; base..., veg_file=ones(length(x)), tstop=60day,
                         refdate="2020-06-01 00:00",
                         veg_water_file=[0.0 0.0 5.0; 4000day 0.0 5.0])
    integ = init(AeolisProblem(w))
    solve!(integ)
    @test integ.u.veg_fW[1, 2] == 0 && integ.u.veg_s[1, 2] < tr.s_w
end

@testset "marram roots and dune erosion" begin
    e = marram_experiment(:root_erosion)
    @test e.figlus_reduction ≈ 1 / 3.5
    @test e.stand_observed.kosmalla - 0.05 <= e.stand_reduction <=
          e.stand_observed.debattisti + 0.05
    @test Aeolis._resisted_depth(1.0, 0.0, 0.0, 5.0, 2.5) == 1.0
    @test Aeolis._resisted_depth(1.0, 0.0, 2.0, 2.5, 2.5) == 2.0
    @test Aeolis._resisted_depth(1.0, 0.0, 0.5, 2.5, 2.5) == 1.5
    dir = joinpath(@__DIR__, "data", "case3_erosion_avalanching")
    p = read_config(joinpath(dir, "aeolis.txt"); veg_model=:marram, output_file=nothing)
    integ = init(AeolisProblem(p))
    for _ in 1:3
        step!(integ)
    end
    s0 = copy(integ.u)
    m = integ.model
    m.dt_opt = 86400.0         # a day of collision, so the cut spans many cells
    function eroded(cr, zr)
        s = copy(s0)
        s.veg_cr .= cr
        s.veg_zr .= zr
        Aeolis.dune_erosion!(s, m, integ.t)
        return sum(s0.zb[:, 2] .- s.zb[:, 2])
    end
    bare = eroded(0.0, 0.0)
    @test bare > 0
    @test eroded(2.0, 0.0) == bare                       # no rooted layer, no resistance
    @test eroded(1.0, 10.0) / bare ≈ 2.5 / 3.5 atol = 0.06
    @test eroded(2.5, 10.0) / bare ≈ 0.5 atol = 0.06
end
