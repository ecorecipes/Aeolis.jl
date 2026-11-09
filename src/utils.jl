# Numerical helpers shared by the processes. Several reproduce a numpy or scipy routine
# exactly where the model's results depend on its conventions (index placement, boundary
# handling, rounding of grid lengths); each says which.

"""
    rotate(x, y, alpha; origin=(0, 0))

Rotate the points `(x, y)` by `alpha` degrees about `origin`, in the sense AeoLiS uses: the
row vector `[x y]` is multiplied by `[cos a -sin a; sin a cos a]`, so a positive angle turns
the points clockwise.
"""
function rotate(x::AbstractArray, y::AbstractArray, alpha::Real; origin=(0.0, 0.0))
    a = alpha / 180.0 * pi
    c, s = cos(a), sin(a)
    xr = x .- origin[1]
    yr = y .- origin[2]
    return (xr .* c .+ yr .* s .+ origin[1], .-xr .* s .+ yr .* c .+ origin[2])
end
function rotate(x::Real, y::Real, alpha::Real; origin=(0.0, 0.0))
    a = alpha / 180.0 * pi
    c, s = cos(a), sin(a)
    xr, yr = x - origin[1], y - origin[2]
    return (xr * c + yr * s + origin[1], -xr * s + yr * c + origin[2])
end

"""
    rotate!(x, y, alpha; origin=(0, 0))

[`rotate`](@ref) in place: `x` and `y` are overwritten by the rotated points, with the
same arithmetic and so the same values.
"""
function rotate!(x::AbstractArray, y::AbstractArray, alpha::Real; origin=(0.0, 0.0))
    a = alpha / 180.0 * pi
    c, s = cos(a), sin(a)
    o1, o2 = origin
    for k in eachindex(x, y)
        xr = x[k] - o1
        yr = y[k] - o2
        x[k] = xr * c + yr * s + o1
        y[k] = -xr * s + yr * c + o2
    end
    return x, y
end

"""
    normalize_along(x, dim; ref=nothing, fill=0.0)

`x` divided by its sum along `dim` (or by `ref`), with `fill` where the divisor is zero
(numpy-style `normalize` of AeoLiS).
"""
function normalize_along(x::AbstractArray, dim::Integer; ref=nothing, fill=0.0)
    r = ref === nothing ? sum(x; dims=dim) : ref
    return broadcast((a, b) -> b != 0 ? a / b : oftype(a, fill), x, r)
end

"""
    prevent_tiny_negatives!(x, max_error=1e-10, replacement=0.0)

Replace the values of `x` in `(-max_error, 0)` by `replacement`.
"""
function prevent_tiny_negatives!(x::AbstractArray, max_error=1e-10, replacement=0.0)
    @inbounds for i in eachindex(x)
        if x[i] < 0 && x[i] > -max_error
            x[i] = replacement
        end
    end
    return x
end

"""
    apply_mask!(arr, mask)

Multiply `arr` by the real part of the complex `mask` and add its imaginary part. A mask
covering the first dimensions of `arr` is repeated over the rest; a scalar mask applies to
every cell.
"""
function apply_mask!(arr::AbstractArray, mask::Number)
    arr .*= real(mask)
    arr .+= imag(mask)
    return arr
end
function apply_mask!(arr::AbstractArray, mask::AbstractArray)
    m = reshape(mask, size(mask)..., ntuple(_ -> 1, ndims(arr) - ndims(mask))...)
    arr .*= real.(m)
    arr .+= imag.(m)
    return arr
end

# ---------------------------------------------------------------------------------------
# Time series

"""
    CircularSeries(t, values)

A linear interpolant of a time series that repeats beyond its range the way AeoLiS's
`interp_circular` does: a time `x` maps to `tmin + mod(x - tmax - 1, tmax - tmin + 1)` and is
then interpolated linearly, clamped at the ends (`numpy.interp`).
"""
struct CircularSeries{I}
    tmin::Float64
    tmax::Float64
    itp::I
end
function CircularSeries(t::AbstractVector, v::AbstractVector)
    itp = DataInterpolations.LinearInterpolation(Float64.(v), Float64.(t);
                                                 extrapolation=DataInterpolations.ExtrapolationType.Constant)
    return CircularSeries(Float64(minimum(t)), Float64(maximum(t)), itp)
end
function (c::CircularSeries)(x::Real)
    xr = c.tmax - c.tmin
    xm = c.tmin + mod(x - c.tmax - 1.0, xr + 1.0)
    return c.itp(xm)
end

"""
    CircularStepSeries(t, values)

The value of a time series at the last time not after `x`, repeating with period
`tmax - tmin` (AeoLiS's `interp_circular_nearest`, used for wind).
"""
struct CircularStepSeries{I}
    tmin::Float64
    tmax::Float64
    itp::I
end
function CircularStepSeries(t::AbstractVector, v::AbstractVector)
    itp = DataInterpolations.ConstantInterpolation(Float64.(v), Float64.(t); dir=:left,
                                                   extrapolation=DataInterpolations.ExtrapolationType.Constant)
    return CircularStepSeries(Float64(minimum(t)), Float64(maximum(t)), itp)
end
function (c::CircularStepSeries)(x::Real)
    xm = c.tmin + mod(x - c.tmin, c.tmax - c.tmin)
    return c.itp(xm)
end

"""
    clamped_interp(x, xp, fp)

Linear interpolation clamped at the ends (`numpy.interp`).
"""
function clamped_interp(x::Real, xp::AbstractVector, fp::AbstractVector)
    x <= xp[1] && return float(fp[1])
    x >= xp[end] && return float(fp[end])
    i = searchsortedlast(xp, x)
    i >= length(xp) && return float(fp[end])
    xp[i + 1] == xp[i] && return float(fp[i + 1])
    return fp[i] + (x - xp[i]) * (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
end

# ---------------------------------------------------------------------------------------
# Sparse matrices

"""
    diags(diagonals, offsets, n)

An `n` by `n` sparse matrix with `diagonals[k]` on diagonal `offsets[k]`, placed as
`scipy.sparse.diags` places them: for an offset `d >= 0` entry `i` goes to `(i, i+d)`, for
`d < 0` to `(i-d, i)`. The AeoLiS solvers rely on this placement, including where a diagonal
is taken from a shifted copy of a coefficient array.
"""
function diags(diagonals, offsets, n::Integer)
    nnz = sum(d -> max(n - abs(d), 0), offsets; init=0)
    I = Vector{Int}(undef, nnz)
    J = Vector{Int}(undef, nnz)
    V = Vector{Float64}(undef, nnz)
    k = 0
    for (v, d) in zip(diagonals, offsets)
        m = n - abs(d)
        length(v) >= m ||
            throw(ArgumentError("diagonal at offset $d needs $m values, got $(length(v))"))
        for i in 1:m
            k += 1
            I[k], J[k] = d >= 0 ? (i, i + d) : (i - d, i)
            V[k] = v[i]
        end
    end
    resize!(I, k)
    resize!(J, k)
    resize!(V, k)
    return SparseArrays.sparse(I, J, V, n, n)
end

# ---------------------------------------------------------------------------------------
# Grids

"""
    np_arange(start, stop, step)

The points of `numpy.arange(start, stop, step)` for floats: `ceil((stop - start) / step)`
points `start + i * step`.
"""
function np_arange(start::Real, stop::Real, step::Real)
    n = max(0, ceil(Int, (stop - start) / step))
    return [start + i * step for i in 0:(n - 1)]
end

"""
    exact_grid(xmin, xmax, ymin, ymax, dx, dy)

The regular grid of AeoLiS's `get_exact_grid`, as `(x, y)` arrays of size `(nx, ny)`.
"""
function exact_grid(xmin, xmax, ymin, ymax, dx, dy)
    xs = np_arange(floor(xmin / dx) * dx, ceil(xmax / dx) * dx, dx)
    ys = np_arange(floor(ymin / dy) * dy, ceil(ymax / dy) * dy, dy)
    return ([xi for xi in xs, _ in ys], [yi for _ in xs, yi in ys])
end

"""
    fftfreq(n, d)

The sample frequencies of `numpy.fft.fftfreq(n, d)`.
"""
function fftfreq(n::Integer, d::Real)
    val = 1.0 / (n * d)
    k = vcat(0:((n - 1) ÷ 2), (-(n ÷ 2)):-1)
    return k .* val
end

"""
    pad_odd_reflect(v, w)

`v` padded by `w` points at each end with `numpy.pad(..., 'reflect', reflect_type='odd')`,
which reflects repeatedly when `w` exceeds `length(v) - 1`.
"""
function pad_odd_reflect(v::AbstractVector, w::Integer)
    n = length(v)
    out = zeros(Float64, n + 2w)
    out[(w + 1):(w + n)] .= v
    left, right = w, w
    while left > 0 || right > 0
        old = length(out) - left - right - 1
        if left > 0
            chunk = min(old, left)
            edge = out[left + 1]
            for k in 1:chunk
                out[left + 1 - k] = 2 * edge - out[left + 1 + k]
            end
            left -= chunk
        end
        if right > 0
            chunk = min(old, right)
            e = length(out) - right
            edge = out[e]
            for k in 1:chunk
                out[e + k] = 2 * edge - out[e - k]
            end
            right -= chunk
        end
    end
    return out
end

"""
    pad_edge(a, w)

The matrix `a` padded by `w` cells on every side with its edge values (`numpy.pad` with
`'edge'`).
"""
function pad_edge(a::AbstractMatrix, w::Integer)
    nx, ny = size(a)
    return [a[clamp(i - w, 1, nx), clamp(j - w, 1, ny)]
            for i in 1:(nx + 2w), j in 1:(ny + 2w)]
end

"""
    regular_grid_interp(xg, yg, z, xq, yq, fill)

Bilinear interpolation of `z[i, j]` given at `(xg[i], yg[j])` to the points `(xq, yq)`, with
`fill` outside the grid, as `scipy.interpolate.RegularGridInterpolator` with
`bounds_error=False` does: the cell is found with a left-sided search, clamped to the last
interval.
"""
function regular_grid_interp(xg::AbstractVector, yg::AbstractVector, z::AbstractMatrix,
                             xq::AbstractArray, yq::AbstractArray, fill::Real)
    nxg, nyg = length(xg), length(yg)
    out = similar(xq, Float64)
    @inbounds for k in eachindex(xq, yq)
        x, y = xq[k], yq[k]
        if !(xg[1] <= x <= xg[end]) || !(yg[1] <= y <= yg[end]) || isnan(x) || isnan(y)
            out[k] = fill
            continue
        end
        i = clamp(searchsortedfirst(xg, x) - 1, 1, nxg - 1)
        j = clamp(searchsortedfirst(yg, y) - 1, 1, nyg - 1)
        tx = (x - xg[i]) / (xg[i + 1] - xg[i])
        ty = (y - yg[j]) / (yg[j + 1] - yg[j])
        out[k] = z[i, j] * (1 - tx) * (1 - ty) + z[i, j + 1] * (1 - tx) * ty +
                 z[i + 1, j] * tx * (1 - ty) + z[i + 1, j + 1] * tx * ty
    end
    return out
end

"""
    linear_fit(x, y)

Slope and intercept of the least-squares line through the points (`numpy.polyfit(x, y, 1)`).
"""
function linear_fit(x::AbstractVector, y::AbstractVector)
    A = hcat(x, ones(length(x)))
    c = A \ y
    return c[1], c[2]
end

"""
    max_diff(arr)

The largest `abs(arr[i] - arr[j]) + abs(i - j)` over pairs (AeoLiS's `maxDiff`).
"""
function max_diff(arr::AbstractVector)
    result = 0.0
    n = length(arr)
    for i in 1:n, j in i:n
        v = abs(arr[i] - arr[j]) + abs(i - j)
        v > result && (result = v)
    end
    return result
end

# ---------------------------------------------------------------------------------------
# Filters

function _gaussian_kernel(sigma::Real, truncate::Real=4.0)
    radius = trunc(Int, truncate * sigma + 0.5)
    x = (-radius):radius
    w = exp.(-0.5 .* (x ./ sigma) .^ 2)
    return w ./ sum(w), radius
end

# scipy.ndimage 'reflect' mode: (d c b a | a b c d | d c b a)
function _reflect_index(i::Int, n::Int)
    n == 1 && return 1
    period = 2n
    k = mod(i - 1, period)
    return k < n ? k + 1 : period - k
end

"""
    gaussian_filter1d(a, sigma, dim)

Gaussian filter of `a` along dimension `dim` with `scipy.ndimage.gaussian_filter1d`
semantics (`mode='reflect'`, `truncate=4`). A zero `sigma` returns a copy.
"""
function gaussian_filter1d(a::AbstractArray, sigma::Real, dim::Integer)
    sigma <= 0 && return copy(a)
    w, r = _gaussian_kernel(sigma)
    out = similar(a, Float64)
    n = size(a, dim)
    for I in CartesianIndices(a)
        acc = 0.0
        i = I[dim]
        for (k, off) in enumerate((-r):r)
            J = Base.setindex(Tuple(I), _reflect_index(i + off, n), dim)
            acc += w[k] * a[J...]
        end
        out[I] = acc
    end
    return out
end

"""
    gaussian_filter(a, sigma)

Gaussian filter of `a` along every dimension (`scipy.ndimage.gaussian_filter`).
"""
function gaussian_filter(a::AbstractArray, sigma::Real)
    sigma <= 0 && return copy(a)
    out = Float64.(a)
    for d in 1:ndims(a)
        out = gaussian_filter1d(out, sigma, d)
    end
    return out
end

"""
    grain_size_percentile(p, mass, percent)

Grain size at `percent` of the cumulative mass of the top layer in every cell (AeoLiS's
`calc_grain_size`), interpolated linearly and extrapolated beyond the coarsest fraction.
"""
function grain_size_percentile(grain_size::AbstractVector, mass::AbstractArray,
                               percent::Real)
    nx, ny = size(mass, 1), size(mass, 2)
    D = zeros(nx, ny)
    diameters = vcat(0.0, grain_size)
    for j in 1:ny, i in 1:nx
        cw = cumsum(vcat(0.0, mass[i, j, 1, :]))
        pct = 100 .* cw ./ maximum(cw)
        D[i, j] = _interp_extrapolate(pct, diameters, percent)
    end
    return D
end

function _interp_extrapolate(xp, fp, x)
    n = length(xp)
    k = searchsortedfirst(xp, x)
    lo = clamp(k - 1, 1, n - 1)
    return fp[lo] + (x - xp[lo]) * (fp[lo + 1] - fp[lo]) / (xp[lo + 1] - xp[lo])
end

"""
    mean_grain_size(grain_size, mass)

Mass-weighted mean grain size of the top layer in every cell.
"""
function mean_grain_size(grain_size::AbstractVector, mass::AbstractArray)
    nx, ny = size(mass, 1), size(mass, 2)
    D = zeros(nx, ny)
    for j in 1:ny, i in 1:nx
        m = mass[i, j, 1, :]
        D[i, j] = sum(grain_size .* m ./ sum(m))
    end
    return D
end
