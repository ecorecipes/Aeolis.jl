using Test
using Aeolis
using Random
using NCDatasets
using NonlinearSolveFirstOrder: NewtonRaphson

@testset "Aeolis" begin
    include("test_parameters.jl")
    include("test_utils.jl")
    include("test_processes.jl")
    include("test_marram.jl")
    include("test_problem.jl")
    include("test_steadystate.jl")
    include("test_reference.jl")
    include("test_netcdf.jl")
    include("test_errors.jl")
    include("test_docstrings.jl")
    include("test_aqua.jl")
end
