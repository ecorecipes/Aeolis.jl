# Vignettes

The vignettes run from simple to advanced and together demonstrate what the package does.

| Vignette | Topic |
|---|---|
| [01 Getting started](01_getting_started/01_getting_started.md) | parameters, problems, solutions and plots on a dry transect; fetch and erosion |
| [02 AeoLiS configuration files](02_configuration_files/02_configuration_files.md) | reading and writing `aeolis.txt`, `run_config`, netCDF output, agreement with aeolis-python |
| [03 Transport formulations](03_transport_formulations/03_transport_formulations.md) | thresholds, six saturated-flux formulations, grain velocities |
| [04 Fetch and supply limitation](04_fetch_and_supply/04_fetch_and_supply.md) | the fetch ODE with OrdinaryDiffEq, numerical diffusion, a non-erodible layer |
| [05 Tides, waves and moisture](05_tides_waves_moisture/05_tides_waves_moisture.md) | water levels, runup, infiltration and the moisture threshold |
| [06 Sediment sorting and armouring](06_sediment_sorting/06_sediment_sorting.md) | several fractions in bed layers, coarsening and the fall of supply |
| [07 Vegetation and foredunes](07_vegetation/07_vegetation.md) | Raupach and Okin shear reduction, trapping, growth and burial |
| [08 Sand fences](08_sand_fences/08_sand_fences.md) | sheltering, deposition and burial of a fence |
| [09 Dunes in two dimensions](09_dunes_in_two_dimensions/09_dunes_in_two_dimensions.md) | spectral shear, the separation bubble, a migrating hill |
| [10 Numerical schemes and solvers](10_schemes_and_solvers/10_schemes_and_solvers.md) | implicit, Crank-Nicolson and explicit schemes; trunk, pieter and steady-state solvers; LinearSolve algorithms |
| [11 Callbacks and ensembles](11_callbacks_and_ensembles/11_callbacks_and_ensembles.md) | a nourishment and a storm with SciMLBase callbacks, tstops, an ensemble sweep |
| [12 Calibration](12_calibration/12_calibration.md) | fitting transport parameters to flux observations with Optimization.jl |

Each vignette lives in its own directory `NN_name/NN_name.qmd` and is rendered with
[Quarto](https://quarto.org) using the native Julia engine to three formats:
HTML (`NN_name.html`, self-contained), GitHub-flavoured Markdown (`NN_name.md` plus
`NN_name_files/`) and PDF (`NN_name.pdf`):

```sh
cd vignettes
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
quarto render            # renders every */*.qmd to HTML, GFM and PDF
quarto render --to gfm   # then GFM again: a multi-format render deletes NN_name_files/
```

All three rendered outputs are committed. `julia scripts/sync_vignettes.jl` copies the
Markdown into `docs/src/tutorials/` for Documenter, and CI runs
`sync_vignettes.jl --check` and fails if the copies are stale.

Every `.qmd` must carry `engine: julia` and list `pdf: default` under `format:` in its
front matter. Make sure the package is precompiled before rendering, otherwise Julia's
precompilation output is captured into the first cell of the rendered file.

## PDF

The PDF format (configured in `_quarto.yml`) uses `lualatex` with STIX Two Text/Math for
body text and mathematics and [JuliaMono](https://juliamono.netlify.app) (SIL OFL) for
code. JuliaMono is loaded by fontspec from `../../fonts/JuliaMono/` relative to this
directory. To render elsewhere, download
[`JuliaMono-ttf.tar.gz`](https://github.com/cormullion/juliamono/releases/latest) and unpack
it there, or point `monofontoptions: Path=` at an absolute directory. TeX requirements:
`lualatex`, `koma-script`, `fvextra`, `stix2-otf`.
