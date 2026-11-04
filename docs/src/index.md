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

## Array layout

Spatial arrays are `(nx+1, ny+1)`, with layers and fractions after them: x is the first
dimension. A transect (a vector `xgrid_file`) becomes a grid of three identical rows, as in
AeoLiS; [`transect`](@ref) returns the middle row of saved output.
