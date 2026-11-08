# The calibration experiments of the marram grass model: published experiments and field
# surveys re-run as single-cell simulations, each returned next to the published numbers.

const MARRAM_EXPERIMENTS = (:reijers2021, :nolet2018, :dixon2004, :aptekar2000,
                            :lammers2024, :homberger2024, :root_erosion)

"""
    marram_experiment(name, tr=MarramTraits())

Re-run a published experiment on marram grass with the single-cell model, and return the
model's results next to the published ones. The experiments are:

- `:reijers2021`: the one-year mesocosm of [Reijers2021](@citet). Four tussocks (32 shoots)
  in 6.25 m² of fresh sand were buried by 0, 2 or 4 cm every fortnight from late June 2016,
  that is 0, 0.26 or 0.52 m per growing season. Compares shoots per mesocosm, the
  root:shoot ratio and shoot elongation.
- `:nolet2018`: the growth response to burial over the 2015 growing season of a foredune
  planted in 2011 ([Nolet2018](@citet)). Compares the model's response of a four-year-old
  stand, normalised to its peak, with the fitted quantile curves of NDVI and cover.
- `:dixon2004`: 28 days without water, then rewatering ([Dixon2004](@citet)). Root-zone
  moisture is prescribed to be below wilting. Compares survival.
- `:aptekar2000`: survival after 1, 3, 7 and 13 days of continuous submergence in seawater,
  against the bud viability of rhizomes relative to the controls ([Aptekar2000](@citet)).
- `:lammers2024`: four weeks of a hot July without rain, under open-sided tents, on
  transplants in an embryo dune ([Lammers2024](@citet)). Compares shoot growth with the
  ambient treatment (12.8 mm of rain) for a sparse planting of 10 shoots m⁻² (an
  assumption); the paper found no difference.
- `:homberger2024`: the 132-day pot experiment of [Homberger2024](@citet), with five marram
  shoots per pot of 0.081 m². Watering from 6.5 to 24.4 mm per week raised mean moisture by
  2.2 % by volume per step. Plants died in pots with a mean of 3.3-4.8 %, and biomass rose
  linearly with moisture.

  The model prescribes root-zone moisture at seven levels from 4 % (an assumption within
  the range where plants died) in steps of 2.2 %. Pot hydrology is not modelled: the pots
  held water above the field capacity of free-draining sand. Compares shoot growth with the
  observed death in the driest pots and the continued rise with moisture.

- `:root_erosion`: the reduction of wave-driven dune erosion by roots, `c_r / (tau0 + c_r)`
  (see `dune_erosion!`). It is computed in two ways:
  - with the root reinforcement of about 1 kPa that [Figlus2014](@citet) measured in shear
    tests of their rooted flume dunes, where erosion fell by 33 % under regular and 8 % under
    irregular waves;
  - with the reinforcement of a model stand grown for five years under 0.3 m of burial per
    season, against the 36 % of intact cores of mature *Ammophila* under swash
    ([DeBattisti2020](@citet)) and up to 31 % in scaled flume experiments with marram
    ([Kosmalla2026](@citet)).

Returns a named tuple with `name`, `description`, and model and observed values whose fields
depend on the experiment.
"""
function marram_experiment(name::Symbol, tr::MarramTraits=MarramTraits())
    name in MARRAM_EXPERIMENTS ||
        throw(ArgumentError("unknown marram experiment :$name; the experiments are " *
                            join((":" * String(n) for n in MARRAM_EXPERIMENTS), ", ")))
    return _experiment(Val(name), tr)
end

function _experiment(::Val{:reijers2021}, tr::MarramTraits)
    area = 6.25
    zdot = (0.0, 0.02 / 14, 0.04 / 14)
    runs = map(zdot) do z
        N0 = 32 / area
        c0 = MarramCell(; N=N0, h=0.569, Br=tr.rs0 * tr.m_shoot * N0, zr=0.3, s=tr.s_fc)
        sol = simulate_marram(c0, (t, d) -> (; dzb=z, s=tr.s_fc), tr; days=365,
                              start=Dates.Date(2016, 6, 24))
        return (; shoots=sol.N[end] * area, ratio=sol.Br[end] / (tr.m_shoot * sol.N[end]),
                elongation=elongation(seasonal_burial(z, tr), tr), sol)
    end
    return (; name=:reijers2021,
            description="one-year mesocosm buried 0, 2 or 4 cm per fortnight",
            burial=[seasonal_burial(z, tr) for z in zdot],
            shoots=[r.shoots for r in runs], shoots_observed=[400.0, 400.0, 130.0],
            ratio=[r.ratio for r in runs], ratio_observed=[0.5, 1.0, 2.5],
            elongation=[r.elongation for r in runs],
            elongation_observed=[0.003, 0.005, 0.009], solutions=[r.sol for r in runs])
end

# Nolet et al. (2018), Table 5.1: (b, c) of the Gaussian quantile regressions at the 50th,
# 90th, 95th and 99th quantiles, for the change in NDVI and in cover.
const NOLET_NDVI = ((0.317, 0.261), (0.318, 0.282), (0.314, 0.322), (0.297, 0.469))
const NOLET_COVER = ((0.268, 0.079), (0.313, 0.176), (0.325, 0.228), (0.353, 0.418))

function _experiment(::Val{:nolet2018}, tr::MarramTraits; age::Real=4.0)
    burial = collect(0.0:0.02:0.8)
    shape(d) = burial_response(d, tr) *
               vigour(pressure_step(0.0, d, age * DAYS_PER_YEAR, tr), tr)
    raw = shape.(burial)
    g(d, (b, c)) = exp(-(d - b)^2 / (2c^2))
    ndvi = [g(d, q) for d in burial, q in NOLET_NDVI]
    cover = [g(d, q) for d in burial, q in NOLET_COVER]
    return (; name=:nolet2018,
            description="growth response to burial of a four-year-old planted stand",
            burial, response=raw ./ maximum(raw), ndvi_observed=ndvi,
            cover_observed=cover, quantiles=[50, 90, 95, 99])
end

function _experiment(::Val{:dixon2004}, tr::MarramTraits)
    t0 = MarramTraits(tr; r=0.0, m0=0.0)
    c0 = MarramCell(; N=100.0, h=0.5, s=tr.s_fc)
    dry(t, d) = (; dzb=0.0, s=t < 28 ? tr.s_h : tr.s_fc)
    sol = simulate_marram(c0, dry, t0; days=56, start=Dates.Date(2002, 1, 1))
    return (; name=:dixon2004, description="28 days without water, then rewatered",
            survival=sol.N[end] / c0.N, survival_observed=0.8, solution=sol)
end

function _experiment(::Val{:aptekar2000}, tr::MarramTraits)
    t0 = MarramTraits(tr; r=0.0, m0=0.0)
    days = [1.0, 3.0, 7.0, 13.0]
    surv = map(days) do D
        c0 = MarramCell(; N=100.0, h=0.5, s=tr.s_fc)
        sol = simulate_marram(c0, (t, d) -> (; dzb=0.0, submerged=t < D), t0; days=D + 5,
                              dt=0.25)
        return sol.N[end] / c0.N
    end
    return (; name=:aptekar2000, description="continuous submergence in seawater",
            days, survival=surv, survival_observed=[54.3, 64.8, 51.2, 8.5] ./ 63.8)
end

function _experiment(::Val{:lammers2024}, tr::MarramTraits; N::Real=10.0, pet::Real=4.0)
    rain_ambient(t) = mod(t, 7) == 0 ? 12.8 / 4 : 0.0
    run(rain) = simulate_marram(MarramCell(; N, h=0.3, s=tr.s_fc),
                                (t, d) -> (; dzb=0.0, rain=rain(t), pet), tr; days=28,
                                start=Dates.Date(2022, 7, 7))
    amb, gh = run(rain_ambient), run(t -> 0.0)
    return (; name=:lammers2024, description="four weeks of a hot July without rain",
            growth_ratio=gh.N[end] / amb.N[end], growth_ratio_observed=1.0,
            survival=min(gh.N[end] / N, 1.0), survival_observed=1.0,
            s_greenhouse=gh.s[end], s_ambient=amb.s[end], ambient=amb, greenhouse=gh)
end

function _experiment(::Val{:homberger2024}, tr::MarramTraits)
    # seven water treatments, 2.2 % moisture apart; the driest at 4 % by volume, within the
    # 3.3-4.8 % of the pots in which plants died
    theta = [0.040 + 0.022k for k in 0:6]
    N0 = 5 / 0.081
    runs = map(theta) do th
        s = th / tr.porosity
        return simulate_marram(MarramCell(; N=N0, h=0.3, s), (t, d) -> (; dzb=0.0, s), tr;
                               days=132, start=Dates.Date(2020, 6, 20))
    end
    return (; name=:homberger2024,
            description="132-day pot experiment at seven levels of soil moisture",
            moisture=theta, shoots=[r.N[end] / N0 for r in runs],
            mortality_moisture_observed=(0.033, 0.048), solutions=runs)
end

function _experiment(::Val{:root_erosion}, tr::MarramTraits)
    z = 0.3 / season_length(tr)
    sol = simulate_marram(MarramCell(; N=tr.N_new * tr.N_max, h=tr.h_new, zr=tr.zr_new,
                                     s=tr.s_fc), (t, d) -> (; dzb=z, s=tr.s_fc), tr;
                          days=5 * 365, start=Dates.Date(2020, 1, 1))
    cr = sol.cr[end]
    reduction(c) = c / (tr.tau0 + c)
    return (; name=:root_erosion, description="reduction of dune erosion by roots",
            figlus_cr=1.0, figlus_reduction=reduction(1.0),
            figlus_observed=(regular=0.33, irregular=0.08), stand_cr=cr,
            stand_reduction=reduction(cr),
            stand_observed=(debattisti=0.36, kosmalla=0.31), solution=sol)
end
