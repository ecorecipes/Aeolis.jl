# Each directory in data/ holds a case from aeolis-python (inputs and aeolis.txt) and
# reference.nc: the model state after every step of the Python run, in float64.
const DATA = joinpath(@__DIR__, "data")

function reference_state(a::AbstractArray)
    nd = ndims(a)
    k = nd - 3
    return permutedims(a, vcat([k + 1, k + 2], collect(k:-1:1), [nd]))
end

@testset "agreement with aeolis-python" begin
    for case in sort(filter(d -> isfile(joinpath(DATA, d, "reference.nc")), readdir(DATA)))
        @testset "$case" begin
            dir = joinpath(DATA, case)
            p = read_config(joinpath(dir, "aeolis.txt"))
            NCDataset(joinpath(dir, "reference.nc")) do ref
                t = ref["t"][:]
                integ = init(AeolisProblem(p); save_vars=Symbol[], save_stats=Symbol[])
                names = [Symbol(k) for k in keys(ref) if !(k in ("t", "dt"))]
                data = Dict(n => reference_state(Array{Float64}(Array(ref[String(n)])))
                            for n in names)
                for k in 1:(length(t) - 1)
                    step!(integ)
                    @test integ.t ≈ t[k + 1]
                    for n in names
                        b = selectdim(data[n], ndims(data[n]), k + 1)
                        a = getfield(integ.u, n)
                        sc = maximum(abs, filter(isfinite, b); init=0.0)
                        @test isequal(isfinite.(a), isfinite.(b))
                        fin = isfinite.(b)
                        @test maximum(abs.(a[fin] .- b[fin]); init=0.0) <= 1e-8 * sc + 1e-10
                    end
                end
            end
        end
    end
end
