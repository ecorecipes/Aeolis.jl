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

The vignettes in [`vignettes/`](vignettes/) run from simple to advanced; see
[vignettes/README.md](vignettes/README.md) for the list and how they are rendered.

## Citation

If you use Aeolis.jl, please cite AeoLiS:

*de Vries, S., Hallin, C., van IJzendoorn, C., van Westen, B., Cohn, N., Strypsteen, G., Skaden, J., Agrawal, N., &
Garcia Alvarez, M. (2023). AeoLiS (Version 3.0.1) [Computer software]. https://github.com/openearth/aeolis-python*

## Licence

Aeolis.jl is a translation of AeoLiS and is distributed under the same licence, the GNU General Public License,
version 3 or later ([LICENSE](LICENSE)). AeoLiS is copyright the AeoLiS Development Team.
