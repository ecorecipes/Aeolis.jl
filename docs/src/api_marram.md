# Marram grass

A model of marram grass (*Calamagrostis arenaria*, formerly *Ammophila arenaria*) in one
grid cell. Its state is a [`MarramCell`](@ref):
- shoot density and canopy height;
- root and rhizome mass, and the thickness of the rooted sand layer;
- root-zone soil moisture;
- soil-biota pressure;
- the current spells of submergence and drought.

Each process is a scalar function of the [`MarramTraits`](@ref), with its source and
evidence grade in its docstring. [`marram_rates`](@ref) collects every multiplier and
mortality, [`marram_step`](@ref) advances a cell, and [`simulate_marram`](@ref) runs one
cell under prescribed burial, water and flooding.

[`marram_experiment`](@ref) re-runs the published experiments to which the traits were
calibrated, next to the published numbers.

On the grid, `veg_model = marram` (with `process_vegetation`) makes every cell a marram cell;
the default `veg_model = legacy` is AeoLiS's vegetation. The grid model is set up by
[`Aeolis.initialize_marram!`](@ref):
- each time step, [`Aeolis.marram_cover!`](@ref) gives the shear reduction its cover;
- every `veg_interval`, [`Aeolis.marram_update!`](@ref) advances the cells and
  [`Aeolis.marram_colonise!`](@ref) colonises bare ones;
- with `process_dune_erosion`, roots resist the cut-back of [`Aeolis.dune_erosion!`](@ref):
  sand in a cell's rooted layer costs `1 + veg_cr / tau0` times as much of the waves' work
  as bare sand.

The grid model's parameters are:

| Parameter | Default | Meaning |
|---|---|---|
| `veg_model` | `legacy` | `marram` for this model |
| `veg_file` | none | initial shoot density, as a fraction of `N_max` |
| `veg_traits_file` | none | TOML file of traits (see [`marram_traits`](@ref)) |
| `veg_water_file` | none | time (s), rainfall and potential evapotranspiration (mm per day); without it, roots never lack water |
| `veg_interval` | 86400 s | interval between updates, in accelerated time |
| `veg_random` | true | random colonisation; false colonises at the expected time |

`gamma_vegshear`, `vegshear_type`, `veg_sigma` and `veg_min_elevation` act as for the legacy
model. The `veg_*` state variables hold the cells and their diagnostics, and can be saved like
any other.

## Traits

Evidence is **M** (measured on marram, or on *A. breviligulata* where the source says so),
**D** (derived from measurements) or **P** (a prior without a direct measurement). The range
bounds the values for sensitivity analysis.

```@example traits
using Aeolis, Markdown
rows = ["| `$(r.name)` | $(r.default) | $(r.unit) | $(r.evidence) | $(r.range[1])-$(r.range[2]) | $(isempty(r.source) ? "-" : r.source) | $(r.description) |"
        for r in marram_info()]
Markdown.parse(join(["| Trait | Default | Unit | Evidence | Range | Source | Description |",
                     "|---|---|---|---|---|---|---|", rows...], "\n"))
```

```@docs
MarramTraits
marram_traits
marram_info
write_marram_traits
marram_ranges
marram_sample
```

## One cell

```@docs
MarramCell
marram_rates
marram_step
simulate_marram
marram_experiment
```

## The grid model

These functions are not exported.

```@docs
Aeolis.initialize_marram!
Aeolis.marram_cover!
Aeolis.marram_update!
Aeolis.marram_colonise!
```

## Processes

These functions are not exported.

```@docs
Aeolis.season_length
Aeolis.season
Aeolis.seasonal_burial
Aeolis.burial_response
Aeolis.elongation
Aeolis.smother_hazard
Aeolis.water_stress
Aeolis.evapotranspiration
Aeolis.drainage
Aeolis.bucket_step
Aeolis.drought_hazard
Aeolis.flood_survival
Aeolis.flood_hazard
Aeolis.pressure_step
Aeolis.pressure_target
Aeolis.vigour
Aeolis.root_shoot_ratio
Aeolis.root_target
Aeolis.reinforcement
Aeolis.standing_fraction
Aeolis.frontal_area
Aeolis.effective_cover
Aeolis.model_datetime
Aeolis.day_of_year
```
