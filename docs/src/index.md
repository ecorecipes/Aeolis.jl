# Aeolis.jl

A Julia port of [AeoLiS](https://github.com/openearth/aeolis-python) [deVries2023](@cite),
the process-based model of aeolian sediment transport with spatially and temporally varying
sediment availability [Hoonhout2016](@cite).

A simulation is an [`AeolisProblem`](@ref) of [`AeolisParameters`](@ref), solved with the
SciML interface:

```julia
using Aeolis
p = read_config("aeolis.txt")          # or AeolisParameters(; xgrid_file=..., ...)
sol = solve(AeolisProblem(p), EulerBackward(Trunk()))
sol[:zb]                                # bed level at every output time
```

The model reproduces aeolis-python 3.1.0rc1 to round-off; behaviours of AeoLiS that look
unintended are kept so that results agree, and documented where they occur.

## Tutorials

The tutorials are the package's vignettes and run from simple to advanced. HTML and PDF versions, and the Quarto
sources, are in [`vignettes/`](https://github.com/ecorecipes/Aeolis.jl/tree/main/vignettes) in the repository.

| Vignette | Topic |
|---|---|
| [01 Getting started](tutorials/01_getting_started.md) | parameters, problems, solutions and plots on a dry transect; fetch and erosion |
| [02 AeoLiS configuration files](tutorials/02_configuration_files.md) | reading and writing `aeolis.txt`, `run_config`, netCDF output, agreement with aeolis-python |
| [03 Transport formulations](tutorials/03_transport_formulations.md) | thresholds, six saturated-flux formulations, grain velocities |
| [04 Fetch and supply limitation](tutorials/04_fetch_and_supply.md) | the fetch ODE with OrdinaryDiffEq, numerical diffusion, a non-erodible layer |
| [05 Tides, waves and moisture](tutorials/05_tides_waves_moisture.md) | water levels, runup, infiltration and the moisture threshold |
| [06 Sediment sorting and armouring](tutorials/06_sediment_sorting.md) | several fractions in bed layers, coarsening and the fall of supply |
| [07 Vegetation and foredunes](tutorials/07_vegetation.md) | Raupach and Okin shear reduction, trapping, growth and burial |
| [08 Sand fences](tutorials/08_sand_fences.md) | sheltering, deposition and burial of a fence |
| [09 Dunes in two dimensions](tutorials/09_dunes_in_two_dimensions.md) | spectral shear, the separation bubble, a migrating hill |
| [10 Numerical schemes and solvers](tutorials/10_schemes_and_solvers.md) | implicit, Crank-Nicolson and explicit schemes; trunk, pieter and steady-state solvers; LinearSolve algorithms |
| [11 Callbacks and ensembles](tutorials/11_callbacks_and_ensembles.md) | a nourishment and a storm with SciMLBase callbacks, tstops, an ensemble sweep |
| [12 Calibration](tutorials/12_calibration.md) | fitting transport parameters to flux observations with Optimization.jl |
| [13 Marram grass](tutorials/13_marram_grass.md) | a model of marram grass in one cell, and the published experiments it is calibrated to |
| [14 A foredune built by marram grass](tutorials/14_marram_foredune.md) | marram grass on a transect: trapping, burial, spread, and planting density |
| [15 Uncertainty in the marram grass model](tutorials/15_marram_uncertainty.md) | Sobol indices of the priors, history matching against the calibration experiments, and an ensemble of foredunes and storms |

## Processes and formulations

| Process | Formulation |
|---|---|
| Shear perturbation by topography | [Weng1991](@cite) (spectral), [Kroy2002](@cite) (`quasi2d`) |
| Grain velocity | [Duran2011](@cite) |
| Equilibrium transport | [Bagnold1937](@cite) and others, [vanRijn2020](@cite) |
| Runup | [Stockdon2006](@cite) |
| Dune erosion | [Palmsten2012](@cite) |
| Shear reduction by vegetation and roughness | [Raupach1992](@cite), [Raupach1993](@cite), [Okin2008](@cite) |
| Vegetation growth | AeoLiS, or [Duran2013](@cite) |
| Marram grass ([Marram grass](api_marram.md)) | [Nolet2018](@cite), [Reijers2021](@cite), [Homberger2026](@cite), [vanderPutten1993](@cite) and others |

## Array layout

Spatial arrays are `(nx+1, ny+1)`, with layers and fractions after them: x is the first
dimension. A transect (a vector `xgrid_file`) becomes a grid of three identical rows, as in
AeoLiS; [`transect`](@ref) returns the middle row of saved output.
