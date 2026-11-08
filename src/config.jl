# AeoLiS configuration files (`aeolis.txt`, AeoLiS `inout.py`): `key = value` lines with `%`
# comments; values are booleans (T/F), numbers, lists separated by spaces, `None`, or names of
# files that are loaded as arrays.

# Spatial files hold one row per grid row (y); they are transposed to the x-first layout.
const SPATIAL_FILES = (:xgrid_file, :ygrid_file, :bed_file, :ne_file, :veg_file,
                       :fence_file,
                       :threshold_file, :supply_file, :bedcomp_file, :wave_mask, :tide_mask,
                       :runup_mask, :threshold_mask, :gw_mask, :vver_mask)
const SERIES_FILES = (:wind_file, :tide_file, :wave_file, :meteo_file, :veg_water_file)
# Files that are kept as paths rather than loaded as arrays.
const PATH_FILES = (:veg_traits_file,)

"""
    load_array(file; spatial=false)

The numbers in a text file as AeoLiS's `numpy.loadtxt` reads them: a vector for a single row
or column, otherwise a matrix (transposed for spatial files, so that x is the first
dimension). Complex numbers are allowed (masks).
"""
function load_array(file::AbstractString; spatial::Bool=false)
    a = try
        DelimitedFiles.readdlm(file, Float64; comments=true, comment_char=('#'))
    catch err
        err isa ArgumentError || rethrow()
        raw = DelimitedFiles.readdlm(file, String; comments=true, comment_char=('#'))
        parse.(ComplexF64, replace.(raw, "j" => "im"))
    end
    (size(a, 1) == 1 || size(a, 2) == 1) && return vec(a)
    return spatial ? permutedims(a) : a
end

function _parse_scalar(val::AbstractString)
    val == "T" && return true
    val == "F" && return false
    occursin(r"^-?\d+$", val) && return parse(Int, val)
    x = tryparse(Float64, val)
    x !== nothing && return x
    val == "None" && return nothing
    return String(val)
end

"""
    parse_value(val; dir=".", key=nothing)

The value of a configuration entry: everything after `%` is a comment; a value with spaces
is a list; `T`/`F` are booleans; `None` and an empty value are `nothing`; the name of an
existing file (relative to `dir`) is loaded with [`load_array`](@ref).
"""
function parse_value(val::AbstractString; dir::AbstractString=".", key=nothing)
    val = strip(first(split(val, '%')))
    isempty(val) && return nothing
    if occursin(r"\s", val)
        return [_parse_scalar(String(x)) for x in split(val)]
    end
    v = _parse_scalar(String(val))
    if v isa String
        path = isabspath(v) ? v : joinpath(dir, v)
        key in PATH_FILES && return isfile(path) ? normpath(path) : v
        if isfile(path)
            return load_array(path; spatial=key in SPATIAL_FILES)
        end
    end
    return v
end

"""
    read_config(file; kwargs...)

The [`AeolisParameters`](@ref) of an AeoLiS configuration file, with files loaded relative to
the file's directory. Keywords override entries of the file. Keys that AeoLiS ignores (and
`nfractions`, which follows from `grain_size`) are skipped with a warning listing them.
`output_file` defaults to the configuration file name with extension `.nc`, and `h` to `z`.
"""
function read_config(file::AbstractString; kwargs...)
    isfile(file) || throw(ConfigurationError("configuration file not found: $file"))
    dir = dirname(abspath(file))
    vals = Dict{Symbol,Any}()
    ignored = Symbol[]
    for line in eachline(file)
        (occursin('=', line) && !startswith(strip(line), "%")) || continue
        parts = split(line, '=')
        key = Symbol(strip(parts[1]))
        v = parse_value(parts[2]; dir, key)
        if haskey(PARAMETER_INDEX, key)
            T = PARAMETER_INDEX[key][2]
            if v === nothing && !(Nothing <: T)
                continue
            end
            if T in (String, Union{Nothing,String}) && v isa AbstractArray
                v = join(string.(v), " ")
            elseif T in (Vector{Symbol},) && v isa AbstractArray
                v = Symbol.(string.(v))
            elseif T === Vector{Symbol} && v !== nothing
                v = [Symbol(string(v))]
            elseif T in (Vector{Float64}, VecOrMat{Float64}) && v isa AbstractVector
                v = Float64.(v)
            elseif T === InputArray && v isa AbstractString
                throw(ConfigurationError("$key names $(repr(v)), which is not a file in $dir"))
            end
            vals[key] = v
        else
            push!(ignored, key)
        end
    end
    isempty(ignored) ||
        @warn "configuration keys not used by AeoLiS were skipped" file keys = ignored
    get(vals, :output_file, nothing) === nothing &&
        (vals[:output_file] = splitext(abspath(file))[1] * ".nc")
    get(vals, :h, nothing) === nothing && (vals[:h] = get(vals, :z, 10.0))
    merge!(vals, Dict{Symbol,Any}(kwargs))
    return _make_parameters(vals)
end

_format_value(v::Bool) = v ? "T" : "F"
_format_value(v::Nothing) = "None"
_format_value(v::Symbol) = String(v)
_format_value(v::AbstractVector) = join(_format_value.(v), " ")
_format_value(v) = string(v)

"""
    write_config(file, p::AeolisParameters)

Write `p` as an AeoLiS configuration file. Input arrays are written to text files next to it
(`bed_file` to `bed.txt`, and so on) in the layout AeoLiS reads.
"""
function write_config(file::AbstractString, p::AeolisParameters)
    dir = dirname(abspath(file))
    w = maximum(length.(String.(PARAMETER_NAMES)))
    open(file, "w") do io
        println(io, "%"^70)
        println(io, "%% ", rpad("AeoLiS model configuration (written by Aeolis.jl)", 64),
                " %%")
        println(io, "%"^70)
        println(io)
        for name in sort(collect(PARAMETER_NAMES))
            v = getfield(p, name)
            if v isa AbstractArray && (name in SPATIAL_FILES || name in SERIES_FILES)
                fname = replace(String(name), "_file" => "") * ".txt"
                a = name in SPATIAL_FILES && ndims(v) == 2 ? permutedims(v) : v
                DelimitedFiles.writedlm(joinpath(dir, fname),
                                        ndims(a) == 1 ? permutedims(a) : a, ' ')
                println(io, rpad(String(name), w), " = ", fname)
            else
                println(io, rpad(String(name), w), " = ", _format_value(v))
            end
        end
    end
    return file
end

"""
    run_config(file; outputfile=nothing, kwargs...)

Run the AeoLiS model of a configuration file, like `aeolis run file`: solve with the scheme
and solver of the file and, if NCDatasets is loaded, write the output to `outputfile` (the
file's `output_file`). Returns the [`AeolisSolution`](@ref). Keywords go to `solve`.
"""
function run_config(file::AbstractString; outputfile=nothing, kwargs...)
    p = read_config(file)
    prob = AeolisProblem(p)
    sol = solve(prob; kwargs...)
    out = outputfile === nothing ? p.output_file : outputfile
    if out !== nothing
        ext = Base.get_extension(@__MODULE__, :AeolisNCDatasetsExt)
        if ext === nothing
            @warn "load NCDatasets to write netCDF output" outputfile = out
        else
            write_netcdf(out, sol)
        end
    end
    return sol
end

"""
    write_netcdf(file, sol::AeolisSolution)

Write the saved output of `sol` to a netCDF file with the variables, dimensions and
attributes of AeoLiS output (transects as one row). Needs NCDatasets.jl to be loaded.
"""
function write_netcdf end
