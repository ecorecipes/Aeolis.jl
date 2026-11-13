"""
    Aeolis

A Julia port of AeoLiS, the process-based model of supply-limited aeolian sediment transport
(de Vries et al. 2023): wind and its perturbation by topography, the shear velocity threshold
with moisture, salt, sheltering and non-erodible layers, equilibrium transport, multiple
sediment fractions in bed layers, tides, waves, runup, groundwater and moisture,
avalanching, vegetation, sand fences and dune erosion.

A model is an [`AeolisProblem`](@ref) built from [`AeolisParameters`](@ref) (or an AeoLiS
configuration file, [`read_config`](@ref)) and solved with the SciML interface (`init`,
`step!`, `solve`) by a scheme ([`EulerBackward`](@ref), [`CrankNicolson`](@ref),
[`EulerForward`](@ref)) and a transport solver ([`Trunk`](@ref), [`Pieter`](@ref),
[`SteadyState`](@ref), [`SteadyStatePieter`](@ref), or [`SciMLSteadyState`](@ref) with a
SciML steady-state solver). Loading NCDatasets enables [`write_netcdf`](@ref).
"""
module Aeolis

using CommonSolve: CommonSolve, init, solve, solve!, step!
using SciMLBase: SciMLBase, CallbackSet, DiscreteCallback, EnsembleProblem, ReturnCode,
                 remake, reinit!, terminate!, u_modified!, add_tstop!
using DataInterpolations: DataInterpolations
using Dates: Dates
using DelimitedFiles: DelimitedFiles
using FFTW: FFTW
using LinearAlgebra: LinearAlgebra
using LinearSolve: LinearSolve
using PrecompileTools: PrecompileTools
using Random: Random
using RecipesBase: RecipesBase
using Serialization: Serialization
using SparseArrays: SparseArrays
using SpecialFunctions: SpecialFunctions
using Statistics: mean
using TOML: TOML

include("errors.jl")
include("parameters.jl")
include("state.jl")
include("utils.jl")
include("grid.jl")
include("bed.jl")
include("wind.jl")
include("shear.jl")
include("rotation.jl")
include("threshold.jl")
include("transport.jl")
include("hydro.jl")
include("avalanching.jl")
include("vegetation.jl")
include("marram.jl")
include("marram_experiments.jl")
include("marram_grid.jl")
include("fences.jl")
include("erosion.jl")
include("solvers.jl")
include("model.jl")
include("problem.jl")
include("config.jl")
include("recipes.jl")
include("precompile.jl")

export AeolisError, ConfigurationError, UnknownOptionError, UnsupportedOptionError
export AeolisParameters, parameter_info, remake_parameters
export AeolisState, state_info
export MarramTraits, marram_traits, marram_info, MarramCell, marram_rates, marram_step,
       simulate_marram, marram_experiment, marram_ranges, marram_sample,
       write_marram_traits
export AeolisProblem, AeolisIntegrator, AeolisSolution, transect, statistic
export AbstractAeolisAlgorithm, EulerBackward, CrankNicolson, EulerForward,
       default_algorithm
export AbstractTransportSolver, Trunk, Pieter, SteadyState, SteadyStatePieter,
       SciMLSteadyState
export read_config, write_config, run_config, write_netcdf, load_array
export save_restart
export init, solve, solve!, step!, remake, reinit!, terminate!, u_modified!, add_tstop!
export CallbackSet, DiscreteCallback, EnsembleProblem, ReturnCode

end
