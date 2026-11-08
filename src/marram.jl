# Marram grass (Calamagrostis arenaria, formerly Ammophila arenaria) as a model of one grid
# cell: shoot density, canopy height, roots and rhizomes, the rooted sand layer, root-zone
# soil water, soil-biota pressure, and the spells of submergence and drought. Each process
# is a scalar function of the traits, `marram_rates` collects them, `marram_step` advances
# a cell by a time step, and `simulate_marram` runs one cell under prescribed forcing.

const DAYS_PER_YEAR = 365.25

# (name, default, unit, evidence, low, high, source, description)
#
# Evidence: :M measured on marram (or A. breviligulata where the source says so); :D derived
# from measurements; :P a prior without a direct measurement. `low` and `high` bound the
# plausible range for sensitivity analysis; `source` is a key of the bibliography.
const MARRAM_TRAIT_TABLE = Tuple{Symbol,Float64,String,Symbol,Float64,Float64,String,
                                 String}[
                                         # growing season
                                         (:season_start, 91.0, "day of year", :M, 60.0,
                                          121.0, "Nolet2018",
                                          "First day of the growing season (1 April)"),
                                         (:season_end, 273.0, "day of year", :M, 213.0,
                                          304.0, "Nolet2018",
                                          "Last day of the growing season (30 September)"),
                                         (:season_ramp, 14.0, "day", :P, 0.0, 30.0,
                                          "Bonte2021",
                                          "Length of the linear ramps at each end of the season"),
                                         # burial
                                         (:burial_window, 30.0, "day", :P, 14.0, 90.0,
                                          "Nolet2018",
                                          "Time scale of the moving average of the rate of bed-level change"),
                                         (:burial_opt, 0.31, "m/season", :M, 0.27, 0.35,
                                          "Nolet2018",
                                          "Burial per growing season at which growth is fastest"),
                                         (:burial_c_lo, 2.0, "m/season", :D, 1.0, 5.0,
                                          "Reijers2021",
                                          "Width of the burial response below the optimum"),
                                         (:burial_c_hi, 0.26, "m/season", :M, 0.08, 0.47,
                                          "Nolet2018",
                                          "Width of the burial response above the optimum"),
                                         # shoots
                                         (:r, 6.3, "1/year", :D, 4.0, 8.0, "Reijers2021",
                                          "Relative growth rate of shoot density per year of full growing season"),
                                         (:N_max, 470.0, "1/m^2", :M, 250.0, 1000.0,
                                          "Kosmalla2025",
                                          "Shoot density of a closed stand"),
                                         (:N_min, 0.5, "1/m^2", :P, 0.1, 2.0, "",
                                          "Shoot density below which a cell is bare"),
                                         (:m0, 0.2, "1/year", :D, 0.1, 0.4, "Kosmalla2025",
                                          "Background shoot mortality"),
                                         (:h_min, 0.1, "m", :P, 0.05, 0.2,
                                          "MaunLapierre1984",
                                          "Exposed canopy height below which buried shoots are smothered"),
                                         (:tau_sm, 56.0, "day", :D, 30.0, 120.0,
                                          "Reijers2021",
                                          "Time scale of death of fully smothered shoots"),
                                         # canopy height
                                         (:h_max, 0.8, "m", :M, 0.65, 1.0, "Kosmalla2025",
                                          "Height of a mature canopy"),
                                         (:e0, 0.003, "m/day", :D, 0.002, 0.004,
                                          "Reijers2021",
                                          "Shoot elongation without burial"),
                                         (:e1, 0.011, "(m/day)/(m/season)", :D, 0.008,
                                          0.014, "Reijers2021",
                                          "Increase of shoot elongation with burial"),
                                         # root-zone water
                                         (:porosity, 0.42, "-", :M, 0.35, 0.45,
                                          "Homberger2026", "Porosity of the root zone"),
                                         (:Z_r, 0.5, "m", :M, 0.3, 1.0, "Homberger2026",
                                          "Depth of the root zone"),
                                         (:s_h, 0.01, "-", :M, 0.005, 0.024,
                                          "Homberger2026",
                                          "Relative soil moisture at the hygroscopic point"),
                                         (:s_w, 0.096, "-", :M, 0.024, 0.114,
                                          "Homberger2024",
                                          "Relative soil moisture at wilting of marram"),
                                         (:s_star, 0.145, "-", :M, 0.12, 0.155,
                                          "Homberger2026",
                                          "Relative soil moisture below which marram closes its stomata"),
                                         (:s_fc, 0.16, "-", :M, 0.08, 0.2, "Homberger2026",
                                          "Relative soil moisture at field capacity"),
                                         (:E_w, 0.1, "mm/day", :M, 0.05, 0.2,
                                          "Homberger2026",
                                          "Evaporation from the root zone at wilting"),
                                         (:E_max, 5.0, "mm/day", :M, 2.3, 19.7,
                                          "Homberger2026",
                                          "Transpiration of a closed canopy without water stress"),
                                         (:K_s, 5040.0, "mm/day", :M, 1000.0, 10000.0,
                                          "Homberger2026",
                                          "Saturated hydraulic conductivity"),
                                         (:b_leak, 4.05, "-", :M, 4.05, 12.1,
                                          "Homberger2026",
                                          "Exponent of the drainage curve"),
                                         (:interception, 0.5, "mm/day", :M, 0.0, 1.0,
                                          "Homberger2026",
                                          "Rainfall intercepted by the canopy"),
                                         (:q_stress, 3.0, "-", :M, 1.0, 4.0,
                                          "Homberger2026",
                                          "Exponent of the water stress curve"),
                                         (:t_dry, 21.0, "day", :D, 14.0, 28.0, "Dixon2004",
                                          "Spell below wilting after which drought kills shoots"),
                                         (:k_dry, 0.03, "1/day", :D, 0.01, 0.06,
                                          "Dixon2004",
                                          "Drought mortality after that spell"),
                                         # flooding
                                         (:D50, 9.1, "day", :D, 1.0, 20.0, "Aptekar2000",
                                          "Continuous submergence in seawater that halves survival"),
                                         (:k_flood, 5.3, "-", :D, 2.0, 8.0, "Aptekar2000",
                                          "Steepness of survival against submergence"),
                                         # soil-biota pressure
                                         (:d_f, 0.35, "m/season", :D, 0.3, 0.6, "Nolet2018",
                                          "Burial per growing season that fully relieves soil-biota pressure"),
                                         (:tau_P, 0.9, "year", :D, 0.5, 2.0, "Nolet2018",
                                          "Time scale of the build-up of soil-biota pressure"),
                                         (:P_seed, 0.02, "-", :P, 0.005, 0.1,
                                          "vanderStoel2002",
                                          "Pressure that seeds the build-up in fresh sand"),
                                         (:tau_relief, 0.25, "year", :P, 0.08, 1.0,
                                          "deRooijvanderGoes1996",
                                          "Time scale of relief from soil-biota pressure by fresh sand"),
                                         (:p_max, 0.7, "-", :P, 0.4, 0.9,
                                          "vanderPutten1993",
                                          "Loss of vigour under full soil-biota pressure"),
                                         # lateral spread and establishment (grid model)
                                         (:v_front, 0.4, "m/season", :M, 0.3, 0.5,
                                          "Keijsers2015",
                                          "Speed of the front of a spreading patch"),
                                         (:f_jump, 0.1, "-", :P, 0.0, 0.3, "Bonte2021",
                                          "Fraction of colonisations that jump beyond the neighbours"),
                                         (:jump_distance, 0.6, "m", :P, 0.4, 1.0,
                                          "Bonte2021",
                                          "Distance of a jump (three cells of 0.2 m in Bonte et al.)"),
                                         (:N_new, 0.05, "-", :P, 0.01, 0.1, "",
                                          "Shoot density of a newly colonised cell, relative to N_max"),
                                         (:h_new, 0.2, "m", :P, 0.1, 0.4, "",
                                          "Canopy height of a newly colonised cell"),
                                         (:zr_new, 0.3, "m", :P, 0.1, 0.5, "",
                                          "Rooted layer of a newly colonised cell"),
                                         (:establishment, 0.0, "1/year", :P, 0.0, 0.05,
                                          "Keijsers2016",
                                          "Establishment away from existing plants, per cell"),
                                         # roots, rooted layer and reinforcement
                                         (:tau_r_grow, 0.02, "year", :D, 0.005, 0.05,
                                          "Reijers2021",
                                          "Response time of roots and rhizomes growing toward their target"),
                                         (:tau_r_decay, 2.0, "year", :P, 1.0, 4.0, "",
                                          "Response time of roots and rhizomes above their target (die-back)"),
                                         (:m_shoot, 0.001, "kg", :D, 0.0005, 0.002,
                                          "Reijers2021",
                                          "Dry mass of a shoot"),
                                         (:Br_max, 1.2, "kg/m^2", :P, 0.5, 2.0, "",
                                          "Largest root and rhizome mass"),
                                         (:rs0, 0.5, "-", :M, 0.3, 0.7, "Reijers2021",
                                          "Root:shoot ratio without burial"),
                                         (:rs1, 1.0, "-", :M, 0.7, 1.3, "Reijers2021",
                                          "Root:shoot ratio at burial rs_d1"),
                                         (:rs2, 2.5, "-", :M, 1.8, 3.2, "Reijers2021",
                                          "Root:shoot ratio at burial rs_d2 and above"),
                                         (:rs_d1, 0.26, "m/season", :M, 0.2, 0.26,
                                          "Reijers2021",
                                          "Burial of the intermediate root:shoot ratio"),
                                         (:rs_d2, 0.52, "m/season", :M, 0.3, 0.52,
                                          "Reijers2021",
                                          "Burial of the largest root:shoot ratio"),
                                         (:zr_max, 1.5, "m", :M, 0.6, 2.0, "Homberger2024",
                                          "Deepest rooted layer"),
                                         (:f_fine, 0.5, "-", :P, 0.3, 0.7, "",
                                          "Fraction of root mass in fine roots"),
                                         (:k_wu, 1.2, "-", :M, 1.0, 1.2, "Wu1979",
                                          "Wu-Waldron orientation factor"),
                                         (:k_fb, 0.4, "-", :D, 0.2, 0.6, "Schwarz2013",
                                          "Fibre-bundle correction of the Wu-Waldron model"),
                                         (:T_r, 15.0, "MPa", :P, 10.0, 40.0, "",
                                          "Tensile strength of roots"),
                                         (:rho_root, 350.0, "kg/m^3", :M, 250.0, 450.0,
                                          "WalkerZinnert2022",
                                          "Root tissue density"),
                                         (:zr_min_cr, 0.3, "m", :P, 0.2, 0.5, "",
                                          "Least depth over which roots are spread for reinforcement"),
                                         (:c_max, 5.0, "kPa", :M, 2.0, 6.0, "Figlus2014",
                                          "Largest root reinforcement"),
                                         (:tau0, 2.5, "kPa", :M, 2.0, 3.0, "Erikson2007",
                                          "Shear strength of moist bare sand"),
                                         # drag
                                         (:a_f, 0.026, "m", :M, 0.013, 0.026,
                                          "LaporteFauret2024",
                                          "Frontal width of a shoot per unit height"),
                                         (:kappa, 0.5, "-", :P, 0.3, 1.0, "",
                                          "Share of frontal area not hidden by overlap"),
                                         (:m_beta, 50.0, "-", :M, 46.0, 76.0, "Walter2012",
                                          "Product of Raupach's m and beta for live grass")]

const MARRAM_TRAIT_NAMES = Tuple(first.(MARRAM_TRAIT_TABLE))
const MARRAM_TRAIT_INDEX = Dict(row[1] => row for row in MARRAM_TRAIT_TABLE)

let fields = [:($(row[1])::Float64) for row in MARRAM_TRAIT_TABLE]
    @eval begin
        """
            MarramTraits(; kwargs...)

        The traits of marram grass (*Calamagrostis arenaria*): growing season, burial
        response, shoot growth and mortality, canopy height, root-zone water, flooding,
        soil-biota pressure, lateral spread, roots and reinforcement, and drag. Every
        keyword defaults to the value in [`marram_info`](@ref), which also gives its unit,
        evidence grade (`:M` measured on marram, `:D` derived from measurements, `:P` a
        prior), plausible range and source. [`marram_traits`](@ref) reads traits from a
        TOML file.

        An unknown keyword, or thresholds of soil moisture or season that are out of
        order, raise an `ArgumentError`.

        # Examples

        ```julia
        tr = MarramTraits(; burial_opt=0.35)
        tr.N_max        # 470.0
        ```
        """
        struct MarramTraits
            $(fields...)
        end
    end
end

function _closest_trait_names(name::Symbol, n::Int=3)
    s = lowercase(String(name))
    score(k) = _edit_distance(s, lowercase(String(k)))
    return sort(collect(MARRAM_TRAIT_NAMES); by=score)[1:n]
end

function _unknown_trait(name)
    return ArgumentError("unknown marram trait :$name; closest names are " *
                         join((":" * String(c) for c in _closest_trait_names(name)), ", "))
end

function MarramTraits(; kwargs...)
    for k in keys(kwargs)
        haskey(MARRAM_TRAIT_INDEX, k) || throw(_unknown_trait(k))
    end
    vals = map(MARRAM_TRAIT_TABLE) do row
        v = get(kwargs, row[1], row[2])
        v isa Real ||
            throw(ArgumentError("marram trait :$(row[1]) expects a number, got $(repr(v))"))
        return Float64(v)
    end
    tr = MarramTraits(vals...)
    _check_traits(tr)
    return tr
end

"""
    MarramTraits(tr::MarramTraits; kwargs...)

A copy of the traits `tr` with the given traits replaced.
"""
function MarramTraits(tr::MarramTraits; kwargs...)
    return MarramTraits(; (n => getfield(tr, n) for n in MARRAM_TRAIT_NAMES)..., kwargs...)
end

function _check_traits(tr::MarramTraits)
    0 < tr.s_h < tr.s_w < tr.s_star < tr.s_fc < 1 ||
        throw(ArgumentError("marram traits need 0 < s_h < s_w < s_star < s_fc < 1, got " *
                            "$(tr.s_h), $(tr.s_w), $(tr.s_star), $(tr.s_fc)"))
    1 <= tr.season_start < tr.season_end <= 366 ||
        throw(ArgumentError("marram traits need 1 <= season_start < season_end <= 366"))
    2 * tr.season_ramp <= tr.season_end - tr.season_start + 1 ||
        throw(ArgumentError("marram trait :season_ramp is longer than half the season"))
    0 <= tr.rs_d1 < tr.rs_d2 ||
        throw(ArgumentError("marram traits need 0 <= rs_d1 < rs_d2"))
    for n in (:burial_c_lo, :burial_c_hi, :burial_window, :N_max, :h_max, :porosity, :Z_r,
              :D50, :tau_P, :P_seed, :tau_relief, :tau_r_grow, :tau_r_decay, :tau_sm, :d_f,
              :rho_root, :m_shoot)
        getfield(tr, n) > 0 || throw(ArgumentError("marram trait :$n must be positive"))
    end
    return tr
end

"""
    marram_traits(file)
    marram_traits(p::AeolisParameters)

Marram traits from a TOML file of `name = value` entries, optionally grouped in tables
(`[burial]`, `[water]`, ...; the table names are ignored). Traits not in the file keep their
defaults. For parameters, the file is `p.veg_traits_file`, and the defaults are used when it
is `nothing`.
"""
function marram_traits(file::AbstractString)
    isfile(file) || throw(ConfigurationError("marram traits file not found: $file"))
    d = TOML.parsefile(file)
    kw = Dict{Symbol,Any}()
    for (k, v) in d
        entries = v isa AbstractDict ? v : Dict(k => v)
        for (k2, v2) in entries
            key = Symbol(k2)
            haskey(kw, key) &&
                throw(ConfigurationError("marram trait :$key is given twice in $file"))
            kw[key] = v2
        end
    end
    return MarramTraits(; kw...)
end
function marram_traits(p::AeolisParameters)
    return p.veg_traits_file === nothing ? MarramTraits() : marram_traits(p.veg_traits_file)
end

"""
    marram_info([name])

A table of every marram trait as `(name, default, unit, evidence, range, source,
description)` named tuples, or the entry of one trait. Evidence is `:M` (measured on
marram, or *A. breviligulata* where the source says so), `:D` (derived from
measurements) or `:P` (a prior without a direct measurement); `source` is a key of the
package bibliography.
"""
function marram_info()
    return [(; name=r[1], default=r[2], unit=r[3], evidence=r[4], range=(r[5], r[6]),
             source=r[7], description=r[8]) for r in MARRAM_TRAIT_TABLE]
end
function marram_info(name::Symbol)
    haskey(MARRAM_TRAIT_INDEX, name) || throw(_unknown_trait(name))
    r = MARRAM_TRAIT_INDEX[name]
    return (; name=r[1], default=r[2], unit=r[3], evidence=r[4], range=(r[5], r[6]),
            source=r[7], description=r[8])
end

function Base.show(io::IO, ::MIME"text/plain", tr::MarramTraits)
    changed = [n for n in MARRAM_TRAIT_NAMES if getfield(tr, n) != MARRAM_TRAIT_INDEX[n][2]]
    print(io, "MarramTraits (", length(MARRAM_TRAIT_NAMES), " traits")
    isempty(changed) && return print(io, ", all defaults)")
    print(io, "; changed:")
    for n in changed
        print(io, "\n  ", n, " = ", getfield(tr, n))
    end
    return print(io, ")")
end
Base.show(io::IO, ::MarramTraits) = print(io, "MarramTraits(...)")
function Base.:(==)(a::MarramTraits, b::MarramTraits)
    return all(n -> getfield(a, n) == getfield(b, n), MARRAM_TRAIT_NAMES)
end
function Base.hash(tr::MarramTraits, h::UInt)
    return hash(Tuple(getfield(tr, n) for n in MARRAM_TRAIT_NAMES), h)
end

# --------------------------------------------------------------------------------------
# Calendar
# --------------------------------------------------------------------------------------

const _REFDATE_FORMATS = (Dates.dateformat"yyyy-mm-dd HH:MM:SS",
                          Dates.dateformat"yyyy-mm-dd HH:MM", Dates.dateformat"yyyy-mm-dd")

"""
    model_datetime(p::AeolisParameters, t)

The date and time of model time `t` (seconds since `refdate`). For vegetation, pass the
accelerated time (`t * accfac`).
"""
function model_datetime(p::AeolisParameters, t::Real)
    s = strip(p.refdate)
    for f in _REFDATE_FORMATS
        d = tryparse(Dates.DateTime, s, f)
        d === nothing && continue
        return d + Dates.Millisecond(round(Int64, 1000 * t))
    end
    return throw(ConfigurationError("refdate $(repr(p.refdate)) is not a date of the form yyyy-mm-dd HH:MM"))
end

"""
    day_of_year(d::Dates.TimeType)

The day of the year of `d` as a real number: 1.0 at the start of 1 January, 1.5 at its noon.
"""
function day_of_year(d::Dates.TimeType)
    dt = Dates.DateTime(d)
    frac = Dates.value(dt - Dates.DateTime(Dates.Date(dt))) / 86_400_000
    return Dates.dayofyear(dt) + frac
end

# --------------------------------------------------------------------------------------
# Process functions. Each is a pure function of scalars and the traits.
# --------------------------------------------------------------------------------------

"""
    season_length(tr)

The length of the growing season in days, `season_end - season_start + 1` (183 for 1 April
to 30 September). Burial rates are converted to m per growing season with it, the unit of
[Nolet2018](@citet).
"""
season_length(tr::MarramTraits) = tr.season_end - tr.season_start + 1

"""
    season(doy, tr)

The growing-season weight φ at day of year `doy` (fractional). It is 0 outside
`season_start`-`season_end`, rises linearly over `season_ramp` days at the start, falls over
`season_ramp` days at the end, and is 1 in between. Growth and lateral spread are multiplied
by it; mortality and burial continue all year.

The window follows the NDVI time course of [Nolet2018](@citet) and the April-August growth of
[Bonte2021](@citet) (evidence **M**); the ramps are a prior (**P**).
"""
function season(doy::Real, tr::MarramTraits)
    a, b, w = tr.season_start, tr.season_end + 1, tr.season_ramp
    (doy < a || doy >= b) && return 0.0
    w == 0 && return 1.0
    return clamp(min(doy - a, b - doy) / w, 0.0, 1.0)
end

"""
    seasonal_burial(zdot, tr)

The burial rate `d` in m per growing season of a rate of bed-level change `zdot` (m per
day): `zdot * season_length(tr)`. Negative for erosion.
"""
seasonal_burial(zdot::Real, tr::MarramTraits) = zdot * season_length(tr)

"""
    burial_response(d, tr)

The growth multiplier `G_b` of burial `d` (m per growing season): a split Gaussian
`exp(-(d - b)^2 / (2c^2))`. The optimum is `b = burial_opt`; the width is `c = burial_c_lo`
below the optimum and `burial_c_hi` above it.

The optimum is that of [Nolet2018](@citet) (evidence **M**, one site and one season). The
widths are fitted to the one-year mesocosm of [Reijers2021](@citet) (evidence **D**). It had
as many shoots without burial as at 0.26 m per season, but a third as many at 0.52 m, and a
symmetric curve cannot do both.
"""
function burial_response(d::Real, tr::MarramTraits)
    c = d < tr.burial_opt ? tr.burial_c_lo : tr.burial_c_hi
    return exp(-(d - tr.burial_opt)^2 / (2 * c^2))
end

"""
    elongation(d, tr)

The shoot elongation rate (m per day) at burial `d` (m per growing season):
`e0 + e1 * max(d, 0)`, a fit to the 0.3, 0.5 and 0.9 cm per day of [Reijers2021](@citet) at 0,
0.26 and 0.52 m per season (evidence **D**).
"""
elongation(d::Real, tr::MarramTraits) = tr.e0 + tr.e1 * max(d, 0.0)

"""
    smother_hazard(h, tr)

The mortality (per day) of shoots whose exposed canopy `h` is below `h_min`:
`max(1 - h / h_min, 0) / tau_sm`, rising from 0 at `h_min` to `1 / tau_sm` when the canopy is
covered. Burial kills by covering the leaves: a single burial of about the canopy height is
lethal, and 0.4 m leaves most plants alive ([MaunLapierre1984](@citet);
[GilbertRipley2010](@citet), *A. breviligulata*).

`tau_sm` is fitted to the mesocosm of [Reijers2021](@citet) (evidence **D**). There, 4 cm of
sand on the plants every fortnight, through a winter without elongation, left a third as
many shoots, which the authors attribute to buried leaves. `h_min` is a prior (**P**).
"""
smother_hazard(h::Real, tr::MarramTraits) = max(1 - h / tr.h_min, 0.0) / tr.tau_sm

"""
    water_stress(s, tr)

The water-stress multiplier `f_W = 1 - ζ` at root-zone relative moisture `s`, with
`ζ = ((s_star - s) / (s_star - s_w))^q_stress` between wilting `s_w` and stomatal closure
`s_star`. `f_W = 1` above `s_star` and 0 below `s_w`. These are the form and values of
[Homberger2026](@citet) (evidence **M**).
"""
function water_stress(s::Real, tr::MarramTraits)
    s >= tr.s_star && return 1.0
    s <= tr.s_w && return 0.0
    return 1 - ((tr.s_star - s) / (tr.s_star - tr.s_w))^tr.q_stress
end

"""
    evapotranspiration(s, pet, cover, tr)

The loss of root-zone water to evapotranspiration (mm per day) at relative moisture `s`.

This is the piecewise function of [Laio2001](@citet), as used by [Homberger2026](@citet):
0 below the hygroscopic point `s_h`, rising to `E_w` at wilting `s_w`, then linearly to the
unstressed rate at `s_star`.

The unstressed rate is `E_w` plus the share `cover` of `min(E_max, pet) - E_w`, where `cover`
(0-1) is the shoot density relative to `N_max`. Transpiration per unit ground area scales with
leaf area, while evaporation from bare sand stays near `E_w`, because a dry surface layer
insulates the root zone. This is why the transplants of [Lammers2024](@citet) were unaffected
when the top 10 cm dried out (evidence **D**).
"""
function evapotranspiration(s::Real, pet::Real, cover::Real, tr::MarramTraits)
    etop = tr.E_w + clamp(cover, 0.0, 1.0) * max(min(tr.E_max, pet) - tr.E_w, 0.0)
    s <= tr.s_h && return 0.0
    s <= tr.s_w && return tr.E_w * (s - tr.s_h) / (tr.s_w - tr.s_h)
    s <= tr.s_star && return tr.E_w + (etop - tr.E_w) * (s - tr.s_w) / (tr.s_star - tr.s_w)
    return etop
end

"""
    drainage(s, tr)

Deep drainage from the root zone (mm per day) above field capacity:
`K_s (exp(b (s - s_fc)) - 1) / (exp(b (1 - s_fc)) - 1)` ([Laio2001](@citet)), with the values
of [Homberger2026](@citet). Their text and table disagree on `b_leak` (12.1 or 4.05).
"""
function drainage(s::Real, tr::MarramTraits)
    s <= tr.s_fc && return 0.0
    b = tr.b_leak
    return tr.K_s * expm1(b * (s - tr.s_fc)) / expm1(b * (1 - tr.s_fc))
end

"""
    bucket_step(s, rain, pet, cover, dt, tr; supply=0.0)

Root-zone relative moisture after `dt` days of rainfall `rain` and potential
evapotranspiration `pet` (mm per day), with canopy share `cover`. The balance is
`n Z_r ds/dt = I - E(s) - L(s)`, where `I` is rainfall less `interception`, plus any capillary
`supply` (mm per day).

Rain is added first, and anything above saturation runs off. Losses are then integrated with
steps that change `s` by at most 0.002, so `s` stays between `s_h` and 1 for any `dt`.
"""
function bucket_step(s::Real, rain::Real, pet::Real, cover::Real, dt::Real,
                     tr::MarramTraits; supply::Real=0.0)
    nz = tr.porosity * tr.Z_r * 1000            # storage in mm
    s = min(Float64(s) + (max(rain - tr.interception, 0.0) + supply) * dt / nz, 1.0)
    left = Float64(dt)
    while left > 0
        loss = evapotranspiration(s, pet, cover, tr) + drainage(s, tr)
        loss <= 0 && break
        h = min(left, 0.002 * nz / loss)
        s = max(s - loss * h / nz, tr.s_h)
        left -= h
    end
    return s
end

"""
    drought_hazard(tdry, tr)

Drought mortality (per day) after a spell of `tdry` days with root-zone moisture below
wilting: `k_dry` once the spell is longer than `t_dry`, and 0 before. In [Dixon2004](@citet),
marram kept its leaf water content for about three weeks without water, and 80 % of plants
recovered after four. This hazard gives 81 % (evidence **D**).
"""
drought_hazard(tdry::Real, tr::MarramTraits) = tdry > tr.t_dry ? tr.k_dry : 0.0

"""
    flood_survival(D, tr)

The share of plants surviving `D` days of continuous submergence in seawater,
`1 / (1 + (D / D50)^k_flood)`. It is fitted to the bud viability of rhizomes of
[Aptekar2000](@citet), relative to their controls: 80 % after 7 days and 13 % after 13
(evidence **D**, rhizome fragments rather than whole plants).
"""
flood_survival(D::Real, tr::MarramTraits) = 1 / (1 + (max(D, 0.0) / tr.D50)^tr.k_flood)

"""
    flood_hazard(D, tr)

The instantaneous mortality (per day) after `D` days of continuous submergence, the hazard of
[`flood_survival`](@ref): `k (D/D50)^k / (D (1 + (D/D50)^k))`.
"""
function flood_hazard(D::Real, tr::MarramTraits)
    D <= 0 && return 0.0
    x = (D / tr.D50)^tr.k_flood
    return tr.k_flood * x / (D * (1 + x))
end

"""
    pressure_step(P, d, dt, tr)

Soil-biota pressure after `dt` days at burial `d` (m per growing season), using the exact
solutions of its equations.

Below its target `P_eq(d)`, pressure builds up logistically from the seed `P_seed`:
`dP/dt = (P + P_seed)(P_eq - P) / tau_P`. Root-feeding nematodes colonise fresh sand within
a month ([vanderStoel2002](@citet)), but their populations take years to build up.

Above its target, fresh sand relieves pressure: `dP/dt = (P_eq - P) / tau_relief`. Regular
burial keeps marram vigorous, and plants decline when sand stops arriving
([deRooijvanderGoes1996](@citet)). `tau_relief`, a few months for roots to grow into new
sand, is a prior (**P**).

The slow start of the build-up reconciles two observations:
- one year on fresh sand without burial cost the mesocosm of [Reijers2021](@citet) only
  4 ± 2 % of its shoots;
- four years after planting, the field response of [Nolet2018](@citet) fell by half without
  burial.

`tau_P` is fitted to both (evidence **D**); `P_seed` is a prior (**P**).
"""
function pressure_step(P::Real, d::Real, dt::Real, tr::MarramTraits)
    Peq = pressure_target(d, tr)
    P = clamp(Float64(P), 0.0, 1.0)
    P > Peq && return _relax(P, Peq, dt, tr.tau_relief * DAYS_PER_YEAR)
    A = Peq + tr.P_seed
    x0 = P + tr.P_seed
    e = exp(-A * dt / (tr.tau_P * DAYS_PER_YEAR))
    return clamp(A * x0 / (x0 + (A - x0) * e) - tr.P_seed, 0.0, 1.0)
end

"""
    pressure_target(d, tr)

The soil-biota pressure toward which a stand relaxes at burial `d` (m per growing season):
`clamp(1 - d / d_f, 0, 1)`. Fresh sand lets marram escape its root-feeding nematodes and
pathogens, and pressure builds where sand stops arriving ([vanderPutten1993](@citet);
[deRooijvanderGoes1996](@citet); [vanderStoel2002](@citet)). `d_f` is fitted to the fall of
the field response of [Nolet2018](@citet) below the optimum (evidence **P**).
"""
pressure_target(d::Real, tr::MarramTraits) = clamp(1 - d / tr.d_f, 0.0, 1.0)

"""
    vigour(P, tr)

The vigour `V = 1 - p_max P` under soil-biota pressure `P` (evidence **P**). It multiplies
both the growth rate and the shoot density that the roots can sustain, `V N_max`. A
degenerating stand therefore thins, rather than only growing more slowly, as marram does on
stable dunes ([vanderPutten1993](@citet)).
"""
vigour(P::Real, tr::MarramTraits) = 1 - tr.p_max * clamp(P, 0.0, 1.0)

"""
    root_shoot_ratio(d, tr)

The root:shoot mass ratio toward which roots grow at burial `d` (m per growing season). It is
linear through `rs0`, `rs1` and `rs2` at 0, `rs_d1` and `rs_d2`, and constant outside. Buried
marram puts more into roots: 0.5, 1.0 and 2.5 at 0, 0.26 and 0.52 m per season in
[Reijers2021](@citet) (evidence **M**).
"""
function root_shoot_ratio(d::Real, tr::MarramTraits)
    d <= 0 && return tr.rs0
    d <= tr.rs_d1 && return tr.rs0 + (tr.rs1 - tr.rs0) * d / tr.rs_d1
    d <= tr.rs_d2 &&
        return tr.rs1 + (tr.rs2 - tr.rs1) * (d - tr.rs_d1) / (tr.rs_d2 - tr.rs_d1)
    return tr.rs2
end

"""
    root_target(N, d, tr)

The root and rhizome mass (kg m⁻²) toward which roots relax: `min(ρ_rs(d) m_shoot N, Br_max)`,
with ρ_rs from [`root_shoot_ratio`](@ref).
"""
function root_target(N::Real, d::Real, tr::MarramTraits)
    return min(root_shoot_ratio(d, tr) * tr.m_shoot * N, tr.Br_max)
end

"""
    reinforcement(Br, zr, tr)

Root reinforcement of the rooted layer (kPa). This is the Wu-Waldron estimate `k_wu T_r RAR`,
reduced by the fibre-bundle factor `k_fb` ([Wu1979](@citet); [Schwarz2013](@citet)), and
capped at `c_max`.

The root area ratio RAR is the volume of fine roots per unit area, spread over the rooted
layer: `f_fine Br / (rho_root max(zr, zr_min_cr))`, with the tissue density of
[WalkerZinnert2022](@citet). Measured reinforcement in dune sand is 0.5-5 kPa
([Figlus2014](@citet)). The tensile strength `T_r` of marram roots has not been measured
(evidence **P**).
"""
function reinforcement(Br::Real, zr::Real, tr::MarramTraits)
    rar = tr.f_fine * max(Br, 0.0) / (tr.rho_root * max(zr, tr.zr_min_cr))
    return min(tr.k_fb * tr.k_wu * tr.T_r * 1000 * rar, tr.c_max)
end

"""
    standing_fraction(u, N)

The share of the blade height left standing in wind speed `u` (m s⁻¹) at shoot density `N`
(m⁻²): `clamp(94.45 - 4.12 u + 0.0195 N, 0, 100) / 100`, Eq. 19 of
[LaporteFauret2024](@citet) (evidence **M**).
"""
standing_fraction(u::Real, N::Real) = clamp(94.45 - 4.12 * u + 0.0195 * N, 0.0, 100.0) / 100

"""
    frontal_area(N, h, u, tr)

The frontal area index λ of the canopy: `kappa N a_f h F(u, N)`. It combines the measured
frontal area of marram shoots ([LaporteFauret2024](@citet)), their reconfiguration in wind
([`standing_fraction`](@ref)), and a share `kappa` of area not hidden by overlapping leaves
(evidence **P**).
"""
function frontal_area(N::Real, h::Real, u::Real, tr::MarramTraits)
    return tr.kappa * max(N, 0.0) * tr.a_f * max(h, 0.0) * standing_fraction(u, N)
end

"""
    effective_cover(lam, tr, gamma_vegshear)

The vegetation cover `rhoveg = min(m_beta lam / gamma_vegshear, 1)`. With it, AeoLiS's shear
reduction `1 / sqrt(1 + gamma_vegshear rhoveg)` equals the Raupach drag partition with the
measured `m beta` of live grass ([Walter2012](@citet)), without its `(1 - m sigma lam)` term
and capped as by [Duran2013](@citet). The legacy model is the case where `rhoveg` is the
cover.
"""
function effective_cover(lam::Real, tr::MarramTraits, gamma::Real)
    return min(tr.m_beta * lam / gamma, 1.0)
end

# --------------------------------------------------------------------------------------
# One cell
# --------------------------------------------------------------------------------------

"""
    MarramCell(; N=0, h=0, Br=0, zr=0.3, s=0.16, P=0, Dsub=0, tdry=0, zdot=0)

The state of marram grass in one cell:
- shoot density `N` (m⁻²);
- canopy height above the bed `h` (m);
- root and rhizome dry mass `Br` (kg m⁻²);
- rooted sand layer `zr` (m);
- root-zone relative moisture `s`;
- soil-biota pressure `P` (0-1);
- current spells of submergence `Dsub` and of moisture below wilting `tdry` (days);
- trailing mean rate of bed-level change `zdot` (m per day).
"""
Base.@kwdef struct MarramCell
    N::Float64 = 0.0
    h::Float64 = 0.0
    Br::Float64 = 0.0
    zr::Float64 = 0.3
    s::Float64 = 0.16
    P::Float64 = 0.0
    Dsub::Float64 = 0.0
    tdry::Float64 = 0.0
    zdot::Float64 = 0.0
end

"""
    marram_rates(cell, forcing, tr)

Every multiplier and rate of the marram model for `cell` (a [`MarramCell`](@ref)) under
`forcing`. The forcing is a named tuple with the day of year `doy`, the rate of bed-level
change `dzb` (m per day), and optionally the wind speed `u` (m s⁻¹).

Returns a named tuple:
- `phi`: the season weight ([`season`](@ref));
- `d`: the burial rate from the trailing mean `cell.zdot` (m per growing season);
- `Gb`, `fW`, `V`: the burial, water and vigour multipliers;
- `growth`: the potential relative growth rate `r phi fW Gb V` (per day), with the shoot
  density approaching `V N_max`;
- `m0`, `m_sm`, `m_dr`, `m_in`: background, smothering, drought and flooding mortality
  (per day), and `mortality`, their sum;
- `elong`: shoot elongation (m per day);
- `P_eq`, `Br_eq`: the targets of pressure and roots;
- `cr`, `lam`: root reinforcement (kPa) and frontal area index.
"""
function marram_rates(c::MarramCell, f, tr::MarramTraits)
    phi = season(f.doy, tr)
    d = seasonal_burial(c.zdot, tr)
    Gb = burial_response(d, tr)
    fW = water_stress(c.s, tr)
    V = vigour(c.P, tr)
    growth = tr.r / DAYS_PER_YEAR * phi * fW * Gb * V
    m0 = tr.m0 / DAYS_PER_YEAR
    m_sm = c.N > 0 ? smother_hazard(c.h, tr) : 0.0
    m_dr = drought_hazard(c.tdry, tr)
    m_in = flood_hazard(c.Dsub, tr)
    u = hasproperty(f, :u) ? f.u : 0.0
    return (; phi, d, Gb, fW, V, growth, m0, m_sm, m_dr, m_in,
            mortality=m0 + m_sm + m_dr + m_in, elong=phi * elongation(d, tr),
            P_eq=pressure_target(d, tr), Br_eq=root_target(c.N, d, tr),
            cr=reinforcement(c.Br, c.zr, tr), lam=frontal_area(c.N, c.h, u, tr))
end

# Exact solution of dN/dt = g N (1 - N/K) - m N over t.
function _logistic(N::Float64, g::Float64, m::Float64, K::Float64, t::Float64)
    N <= 0 && return 0.0
    a = g - m
    b = g / K
    abs(a * t) < 1e-12 && return N / (1 + b * N * t)
    e = exp(a * t)
    return a * N * e / (a + b * N * (e - 1))
end

_relax(x, target, t, tau) = target + (x - target) * exp(-t / tau)

"""
    marram_step(cell, forcing, tr, dt)

The [`MarramCell`](@ref) after `dt` days under `forcing`, a named tuple with:
- `doy`, the day of year;
- `dzb`, the rate of bed-level change (m per day);
- optionally `rain` and `pet` (mm per day), which update the root-zone bucket;
- optionally `s`, a prescribed root-zone moisture, used instead of the bucket;
- optionally `submerged`, `true` while the bed is under the sea, or the fraction of the step
  it was under water. Submergence accumulates while it lasts and resets after a step without
  any; repeated tidal flooding therefore counts as exposure to seawater (evidence **P**);
- optionally `u`, the wind speed (m s⁻¹);
- optionally `supply`, a capillary supply from groundwater (mm per day).

Without `rain` and `s`, moisture keeps its value.

The order of the update:
1. The trailing mean of bed-level change.
2. The bucket, then the spells of drought and submergence.
3. The shoot density. Its exact logistic solution has the mortality of the step, where the
   flooding and drought losses are the exact integrals of their hazards over the step.
4. Canopy height, roots (growing quickly toward their target and dying back slowly above it),
   the rooted layer and pressure ([`pressure_step`](@ref); in a bare cell, pressure fades with
   `tau_relief`).

Every relaxation uses its exact solution, so states stay non-negative and bounded for any
`dt`. A cell whose rooted layer is eroded away, or whose density falls below `N_min`, becomes
bare.
"""
function marram_step(c::MarramCell, f, tr::MarramTraits, dt::Real)
    dt = Float64(dt)
    dzb = Float64(f.dzb)
    zdot = _relax(c.zdot, dzb, dt, tr.burial_window)
    cover = c.N / tr.N_max
    s = if hasproperty(f, :s) && !isnan(f.s)
        Float64(f.s)
    elseif hasproperty(f, :rain)
        supply = hasproperty(f, :supply) ? f.supply : 0.0
        bucket_step(c.s, f.rain, hasproperty(f, :pet) ? f.pet : 0.0, cover, dt, tr; supply)
    else
        c.s
    end
    tdry = s < tr.s_w ? c.tdry + dt : 0.0
    wet = hasproperty(f, :submerged) ? clamp(Float64(f.submerged), 0.0, 1.0) : 0.0
    sub = wet > 0
    Dsub = sub ? c.Dsub + wet * dt : 0.0
    mid = MarramCell(c.N, c.h, c.Br, c.zr, s, c.P, c.Dsub, c.tdry, zdot)
    r = marram_rates(mid, f, tr)
    # exact integrals of the drought and flooding hazards over the step
    loss_dr = tr.k_dry * max(tdry - max(c.tdry, tr.t_dry), 0.0)
    loss_in = sub ? log(flood_survival(c.Dsub, tr) / flood_survival(Dsub, tr)) : 0.0
    m = r.m0 + r.m_sm + (loss_dr + loss_in) / dt
    N = _logistic(c.N, r.growth, m, max(r.V, eps()) * tr.N_max, dt)
    # canopy height: elongation toward h_max, lowered directly by burial
    bur = max(dzb, 0.0)
    h = if N <= 0
        0.0
    elseif r.elong > 0
        hstar = tr.h_max * max(1 - bur / r.elong, 0.0)
        max(_relax(c.h, hstar, dt, tr.h_max / r.elong), 0.0)
    else
        max(c.h - bur * dt, 0.0)
    end
    tau_r = c.Br < r.Br_eq ? tr.tau_r_grow : tr.tau_r_decay
    Br = _relax(c.Br, r.Br_eq, dt, tau_r * DAYS_PER_YEAR)
    # without a host, root-feeding soil biota fade as they are relieved by fresh sand
    P = c.N > 0 ? pressure_step(c.P, r.d, dt, tr) :
        _relax(c.P, 0.0, dt, tr.tau_relief * DAYS_PER_YEAR)
    zr = c.N > 0 ? clamp(c.zr + dzb * dt, 0.0, tr.zr_max) : c.zr
    if N < tr.N_min || (c.N > 0 && zr <= 0)
        N, h = 0.0, 0.0
        zr <= 0 && (Br = 0.0)
    end
    return MarramCell(N, h, Br, zr, s, P, Dsub, tdry, zdot)
end

"""
    simulate_marram(cell, forcing, tr; days, start=Dates.Date(2020, 1, 1), dt=1.0)

Run one cell for `days` days from the date `start`, with steps of `dt` days.
- `forcing(t, date)` returns the forcing of [`marram_step`](@ref) at day `t` (from 0) and
  date `date`, without `doy`; the day of year at the middle of the step is added.
- `cell` is a [`MarramCell`](@ref) and `tr` the [`MarramTraits`](@ref).

Returns a named tuple of vectors: the time `t` (days), the `date`, every state variable, and
the derived quantities and multipliers of [`marram_rates`](@ref) (`phi`, `d`, `Gb`, `fW`,
`V`, `growth`, `mortality`, `cr` and `lam`). The first entry is the initial state.

# Examples

```julia
using Dates
tr = MarramTraits()
sol = simulate_marram(MarramCell(; N=5.0, h=0.5), (t, date) -> (; dzb=1.4e-3), tr;
                      days=365, start=Date(2016, 6, 24))
sol.N[end]
```
"""
function simulate_marram(cell::MarramCell, forcing, tr::MarramTraits; days::Real,
                         start::Dates.TimeType=Dates.Date(2020, 1, 1), dt::Real=1.0)
    days > 0 || throw(ArgumentError("days must be positive, got $days"))
    dt > 0 || throw(ArgumentError("dt must be positive, got $dt"))
    n = ceil(Int, days / dt)
    t0 = Dates.DateTime(start)
    ms(x) = Dates.Millisecond(round(Int64, x * 86_400_000))
    keys_state = fieldnames(MarramCell)
    keys_rate = (:phi, :d, :Gb, :fW, :V, :growth, :mortality, :cr, :lam)
    out = Dict{Symbol,Vector{Float64}}(k => Vector{Float64}(undef, n + 1)
                                       for k in (:t, keys_state..., keys_rate...))
    dates = Vector{Dates.DateTime}(undef, n + 1)
    function record!(i, t, c, f)
        out[:t][i] = t
        dates[i] = t0 + ms(t)
        for k in keys_state
            out[k][i] = getfield(c, k)
        end
        r = marram_rates(c, f, tr)
        for k in keys_rate
            out[k][i] = getfield(r, k)
        end
    end
    c = cell
    record!(1, 0.0, c, merge(forcing(0.0, t0), (; doy=day_of_year(t0))))
    for i in 1:n
        t = (i - 1) * dt
        f = merge(forcing(t, t0 + ms(t)), (; doy=day_of_year(t0 + ms(t + dt / 2))))
        c = marram_step(c, f, tr, dt)
        record!(i + 1, i * dt, c, f)
    end
    return (; t=out[:t], date=dates, (k => out[k] for k in keys_state)...,
            (k => out[k] for k in keys_rate)...)
end
