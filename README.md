# Aeolis.jl

[![Build Status](https://github.com/ecorecipes/Aeolis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ecorecipes/Aeolis.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/Aeolis.jl/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

A Julia port of [AeoLiS](https://github.com/openearth/aeolis-python), the process-based model of supply-limited
aeolian sediment transport in coastal environments, with the SciML problem, integrator and solution interface.

AeoLiS was developed at Delft University of Technology and Deltares by the AeoLiS development team
(de Vries et al. 2023; Hoonhout and de Vries 2016). Aeolis.jl follows aeolis-python 3.1.0rc1 process by process,
and reproduces its results to round-off on the AeoLiS regression cases and on cases covering every process.

## Features

- Wind forcing, roughness lengths, and the shear perturbation by topography: the spectral method of Weng et al.
  (1991) with a separation bubble, and the `quasi2d` (Kroy et al. 2002) and `duna2d` methods on wind-aligned grids.
- The shear velocity threshold from grain size, moisture (ten formulations), salt, sheltering by roughness elements
  and a non-erodible layer.
- Equilibrium transport after Bagnold, Kawamura, Lettau, DK, Sauermann and van Rijn and Strypsteen, with grain
  velocities after Duran et al. (2011).
- Multiple sediment fractions in bed layers, with sorting and armouring.
- Tides, waves, runup (Stockdon et al. 2006), groundwater (Boussinesq), surface moisture with soil water retention
  and hysteresis, and salt.
- Avalanching, vegetation (Raupach and Okin shear reduction, growth, burial, germination), sand fences, and dune
  erosion by waves (Palmsten and Holman 2012).
- Beyond AeoLiS, a model of marram grass in one cell, calibrated to published experiments: shoot density, canopy
  height, roots, rooted depth, root-zone water and soil-biota pressure, with the source and evidence grade of every
  trait (`MarramTraits`, `simulate_marram`, `marram_experiment`), and coupled to the grid with `veg_model = marram`:
  the canopy reduces wind shear, burial, flooding and water drive the plants, bare cells are colonised, and roots
  slow wave erosion of the dune.
- Transport solvers: the implicit `Trunk` and finite-volume `Pieter` solvers with Euler backward, Crank-Nicolson
  or Euler forward schemes, and a `SteadyState` sweep; sparse systems solved with LinearSolve.jl.
- `AeolisProblem`, `init`, `step!`, `solve`, `remake`, SciMLBase callbacks and ensembles; AeoLiS configuration
  files (`read_config`, `write_config`, `run_config`); statistics, restart files; netCDF output in the AeoLiS
  layout through NCDatasets.jl; plot recipes.

## Installation

```sh
git clone https://github.com/ecorecipes/Aeolis.jl
cd Aeolis.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Requires Julia ≥ 1.12.

## Example

```julia
using Aeolis

x = collect(0.0:1.0:100.0)
p = AeolisParameters(; xgrid_file=x, bed_file=-1.0 .+ 0.04 .* x,
                     wind_file=[0.0 12.0 270.0; 86400.0 12.0 270.0],
                     tide_file=[0.0 0.0; 86400.0 0.0], process_tide=true,
                     dt=600.0, tstop=6 * 3600.0, output_times=3600.0)
sol = solve(AeolisProblem(p))
transect(sol, :zb)          # bed level along the transect at each output time

sol = run_config("aeolis.txt")   # an AeoLiS configuration file
```

## Vignettes

The vignettes run from simple to advanced:

| Vignette | Topic |
|---|---|
| [01 Getting started](vignettes/01_getting_started/01_getting_started.md) | parameters, problems, solutions and plots on a dry transect; fetch and erosion |
| [02 AeoLiS configuration files](vignettes/02_configuration_files/02_configuration_files.md) | reading and writing `aeolis.txt`, `run_config`, netCDF output, agreement with aeolis-python |
| [03 Transport formulations](vignettes/03_transport_formulations/03_transport_formulations.md) | thresholds, six saturated-flux formulations, grain velocities |
| [04 Fetch and supply limitation](vignettes/04_fetch_and_supply/04_fetch_and_supply.md) | the fetch ODE with OrdinaryDiffEq, numerical diffusion, a non-erodible layer |
| [05 Tides, waves and moisture](vignettes/05_tides_waves_moisture/05_tides_waves_moisture.md) | water levels, runup, infiltration and the moisture threshold |
| [06 Sediment sorting and armouring](vignettes/06_sediment_sorting/06_sediment_sorting.md) | several fractions in bed layers, coarsening and the fall of supply |
| [07 Vegetation and foredunes](vignettes/07_vegetation/07_vegetation.md) | Raupach and Okin shear reduction, trapping, growth and burial |
| [08 Sand fences](vignettes/08_sand_fences/08_sand_fences.md) | sheltering, deposition and burial of a fence |
| [09 Dunes in two dimensions](vignettes/09_dunes_in_two_dimensions/09_dunes_in_two_dimensions.md) | spectral shear, the separation bubble, a migrating hill |
| [10 Numerical schemes and solvers](vignettes/10_schemes_and_solvers/10_schemes_and_solvers.md) | implicit, Crank-Nicolson and explicit schemes; trunk, pieter and steady-state solvers; LinearSolve algorithms |
| [11 Callbacks and ensembles](vignettes/11_callbacks_and_ensembles/11_callbacks_and_ensembles.md) | a nourishment and a storm with SciMLBase callbacks, tstops, an ensemble sweep |
| [12 Calibration](vignettes/12_calibration/12_calibration.md) | fitting transport parameters to flux observations with Optimization.jl |
| [13 Marram grass](vignettes/13_marram_grass/13_marram_grass.md) | a model of marram grass in one cell, and the published experiments it is calibrated to |
| [14 A foredune built by marram grass](vignettes/14_marram_foredune/14_marram_foredune.md) | marram grass on a transect: trapping, burial, spread, and planting density |

Each vignette is also published as a tutorial in the
[documentation](https://ecorecipes.github.io/Aeolis.jl/), and is rendered to HTML, Markdown and PDF in its directory
under [`vignettes/`](vignettes/). The sources are the `.qmd` files there; [vignettes/README.md](vignettes/README.md)
explains how they are rendered.

## Citation

If you use Aeolis.jl, please cite AeoLiS:

*de Vries, S., Hallin, C., van IJzendoorn, C., van Westen, B., Cohn, N., Strypsteen, G., Skaden, J., Agrawal, N., &
Garcia Alvarez, M. (2023). AeoLiS (Version 3.0.1) [Computer software]. https://github.com/openearth/aeolis-python*

## Licence

Aeolis.jl is a translation of AeoLiS and is distributed under the same licence, the GNU General Public License,
version 3 or later ([LICENSE](LICENSE)). AeoLiS is copyright the AeoLiS Development Team.
