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
calibrated, next to the published numbers. The grid model, selected by
`veg_model = marram`, is not implemented yet; the default `veg_model = legacy` is AeoLiS's
vegetation.

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
```

## One cell

```@docs
MarramCell
marram_rates
marram_step
simulate_marram
marram_experiment
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
