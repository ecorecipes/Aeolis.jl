# Spectral wind shear perturbation over topography (AeoLiS `shear.WindShear`), after Weng et
# al. (1991): the bed is interpolated to a regular computational grid aligned with the wind,
# the perturbation is computed in Fourier space with complex modified Bessel functions, a
# separation bubble may smooth the lee sides first, and the stress is interpolated back.

"""
    WindShear(dx, dy, L, l, z0, buffer_width)

The settings of the spectral shear computation: computational cell sizes, the dune length
scale `L`, the inner layer height `l`, the roughness length `z0` (fixed at its initial value,
as in AeoLiS) and the width of the buffer around the rotated grid.

The perturbation of Weng et al. (1991) assumes a uniform roughness. When `z0` is a grid (an
adaptive `method_roughness`), the logarithmic mean of its positive values is used,
`exp(mean(log(z0)))`, since `z0` enters the solution through `log(l / z0)`. AeoLiS has no
working form for this case.
"""
mutable struct WindShear
    dx::Float64
    dy::Float64
    L::Float64
    l::Float64
    z0::Float64
    buffer_width::Float64
    x0::Float64
    y0::Float64
end
WindShear(dx, dy, L, l, z0, bw) = WindShear(dx, dy, L, l, _effective_z0(z0), bw, 0.0, 0.0)

# AeoLiS interpolates a grid of z0 to the computational grid and then indexes it with the
# mask of Fourier modes (`z0new[valid]`), mixing physical and wavenumber space; numpy fails
# to broadcast the result against the inner-layer height, so a spatially varying z0 raises
# a ValueError there. Here the grid is replaced by its logarithmic mean over the cells where
# it is positive and finite (a uniform grid gives its value exactly).
function _effective_z0(z0)
    z0 isa Real && return Float64(z0)
    v = [Float64(z) for z in z0 if isfinite(z) && z > 0]
    isempty(v) &&
        throw(ConfigurationError("the roughness length z0 is nowhere positive; the " *
                                 "spectral shear (method_shear = :fft) needs a positive z0"))
    allequal(v) && return v[1]
    return exp(mean(log, v))
end

# The computational grid (AeoLiS `set_computational_grid`): a box around the input grid,
# aligned with the wind, plus the buffer.
function _computational_grid!(sh, x::AbstractMatrix, y::AbstractMatrix, udir::Real)
    x0, y0 = mean(x), mean(y)
    xc = [x[1, 1], x[1, end], x[end, 1], x[end, end]]
    yc = [y[1, 1], y[1, end], y[end, 1], y[end, end]]
    if any(v -> abs(udir - v) <= 0.001, -1080:90:1079)
        udir -= 0.1
    end
    bW = zeros(4)
    bL = zeros(4)
    mW = mL = 0.0
    for i in 1:4
        mW, bW[i] = linear_fit([xc[i], xc[i] - sin(deg2rad(udir))],
                               [yc[i], yc[i] - cos(deg2rad(udir))])
        mL, bL[i] = linear_fit([xc[i], xc[i] - sin(deg2rad(udir - 90.0))],
                               [yc[i], yc[i] - cos(deg2rad(udir - 90.0))])
    end
    W = abs(max_diff(bW)) / sqrt(mW^2 + 1) + sh.buffer_width * 2
    L = abs(max_diff(bL)) / sqrt(mL^2 + 1) + sh.buffer_width * 2
    sh.x0, sh.y0 = x0, y0
    return exact_grid(x0 - L / 2, x0 + L / 2, y0 - W / 2, y0 + W / 2, sh.dx, sh.dy)
end

# Interpolate z on the grid (x, y) to (xi, yi) (AeoLiS `WindShear.interpolate`): both are
# turned so the grid is axis-aligned, the grid is extended by odd reflection of its
# coordinates and edge values of z (numpy.pad), and RegularGridInterpolator interpolates
# bilinearly with `fill` outside.
function _shear_interpolate(sh, x, y, z, xi, yi, fill)
    dx = x[2, 1] - x[1, 1]
    dy = y[2, 1] - y[1, 1]
    angle = rad2deg(atan(dy / dx))
    (dx <= 0 && dy <= 0) && (angle += 180.0)
    o = (sh.x0, sh.y0)
    x, y = rotate(x, y, angle; origin=o)
    xi, yi = rotate(xi, yi, angle; origin=o)
    if !issorted(y[1, :]) && !issorted(x[:, 1])
        x, y = rotate(x, y, 180.0; origin=o)
        xi, yi = rotate(xi, yi, 180.0; origin=o)
    end
    w = max(size(x)...)
    # first row of the padded x-coordinates (pad along y, then along x)
    r = [pad_odd_reflect(view(x, c, :), w)[1] for c in axes(x, 1)]
    xg = pad_odd_reflect(r, w)
    # first column of the padded y-coordinates
    P0 = Matrix{Float64}(undef, size(y, 1), size(y, 2) + 2w)
    for c in axes(y, 1)
        P0[c, :] .= pad_odd_reflect(view(y, c, :), w)
    end
    yg = [pad_odd_reflect(view(P0, :, k), w)[1] for k in axes(P0, 2)]
    return _padded_interp(xg, yg, z, w, xi, yi, fill)
end

function _padded_interp(xg, yg, z, w, xq, yq, fill)
    nx, ny = size(z)
    nxg, nyg = length(xg), length(yg)
    out = similar(xq, Float64)
    zz(i, j) = z[clamp(i - w, 1, nx), clamp(j - w, 1, ny)]
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
        out[k] = zz(i, j) * (1 - tx) * (1 - ty) + zz(i, j + 1) * (1 - tx) * ty +
                 zz(i + 1, j) * tx * (1 - ty) + zz(i + 1, j + 1) * tx * ty
    end
    return out
end

"""
    wind_shear!(sh, x, y, z, taux, tauy, u0, udir, separation, c, mu_b, taus0, taun0,
                sep_filter_iterations, zsep_y_filter)

The shear stress components over the bed `z` and the height of the separation bubble, on the
input grid (AeoLiS `WindShear.__call__`). Returns `(taux, tauy, hsep)`.
"""
function wind_shear!(sh::WindShear, x, y, z, taux, tauy, u0, udir, separation, c, mu_b,
                     taus0, taun0, sep_filter_iterations, zsep_y_filter)
    u_angle = 270.0 - udir
    xi, yi = _computational_grid!(sh, x, y, udir)
    o = (sh.x0, sh.y0)
    gcx, gcy = rotate(xi, yi, -u_angle; origin=o)
    gcz = _shear_interpolate(sh, x, y, z, gcx, gcy, 0.0)
    gctx = fill(Float64(taus0), size(gcx))
    gcty = fill(Float64(taun0), size(gcx))
    gcx, gcy = rotate(gcx, gcy, u_angle; origin=o)
    gix, giy = rotate(x, y, u_angle; origin=o)
    gctx, gcty = rotate(gctx, gcty, u_angle)
    z_origin = gcz
    if separation
        zsep = _separation(gcx, gcz, sh.dx, sum(gctx), c, mu_b, sep_filter_iterations,
                           zsep_y_filter)
        z_origin = copy(gcz)
        gcz = max.(gcz, zsep)
    end
    dtaux, dtauy = _shear_perturbation(sh, gcz, u0)
    tau = sqrt.(gctx .^ 2 .+ gcty .^ 2)
    for k in eachindex(tau)
        if tau[k] != 0
            tx = gctx[k] / tau[k] + dtaux[k]
            ty = gcty[k] / tau[k] + dtauy[k]
            gctx[k] = tau[k] * tx
            gcty[k] = tau[k] * ty
        end
    end
    gctx .= max.(gctx, 0.0)
    hsepc = nothing
    if separation
        hsepc = gcz .- z_origin
        delta = 1.0 / (0.2 * 0.5)
        f = min.(max.(1.0 .- delta .* hsepc, 0.0), 1.0)
        gctx .*= f
        gcty .*= f
    end
    gcx, gcy = rotate(gcx, gcy, -u_angle; origin=o)
    gix, giy = rotate(gix, giy, -u_angle; origin=o)
    gctx, gcty = rotate(gctx, gcty, -u_angle)
    tx = _shear_interpolate(sh, gcx, gcy, gctx, gix, giy, taus0)
    ty = _shear_interpolate(sh, gcx, gcy, gcty, gix, giy, taun0)
    hsep = separation ? _shear_interpolate(sh, gcx, gcy, hsepc, gix, giy, 0.0) :
           zeros(size(x))
    return tx, ty, hsep
end

# Shear perturbation in Fourier space (AeoLiS `compute_shear`, Weng et al. 1991), with the
# high-frequency filter of `filter_highfrequenies`.
function _shear_perturbation(sh::WindShear, z::AbstractMatrix, u0::Real)
    u0 == 0 && return zeros(size(z)), zeros(size(z))
    nx, ny = size(z)
    kxv = 2pi .* fftfreq(nx, sh.dx)
    kyv = 2pi .* fftfreq(ny, sh.dy)
    kx = [a for a in kxv, _ in kyv]
    ky = [b for _ in kxv, b in kyv]
    hs = FFTW.fft(complex.(z))
    n1, n2 = 1.0, 2.0
    s1 = n1 / log(1.0 / 0.01 - 1.0)
    s2 = -n2 / log(1.0 / 0.99 - 1.0)
    for k in eachindex(hs)
        px = 2pi / sh.dx / abs(kx[k] == 0 ? 1.0 : kx[k])
        py = 2pi / sh.dy / abs(ky[k] == 0 ? 1.0 : ky[k])
        f1 = 1.0 / (1.0 + exp(-(px + n1 - n2) / s1))
        f2 = 1.0 / (1.0 + exp(-(py + n1 - n2) / s2))
        hs[k] *= f1 * f2
    end
    z0 = sh.z0
    L = sh.L / 4.0
    l = sh.l
    for _ in 1:5
        l = 2 * 0.41^2 * L / log(l / z0)
    end
    hm = 1.0
    for _ in 1:5
        hm = L / sqrt(log(hm / z0))
    end
    ul = log(l / z0) / log(hm / z0)
    dtx = zeros(ComplexF64, size(hs))
    dty = zeros(ComplexF64, size(hs))
    for k in eachindex(hs)
        kk = sqrt(kx[k]^2 + ky[k]^2)
        (kk != 0 && kx[k] != 0) || continue
        sigma = sqrt(1im * L * kx[k] * z0 / l)
        kv0 = _besselk(0, 2 * sigma)
        kv1 = _besselk(1, 2 * sigma)
        term = -1 + (2 * log(l / z0) + kk^2 / kx[k]^2) * sigma * (kv1 / kv0)
        dtx[k] = hs[k] * kx[k]^2 / kk * 2 / ul^2 * term
        kv1y = _besselk(1, 2 * sqrt(2.0) * sigma)
        dty[k] = hs[k] * (kx[k] * ky[k]) / kk * 2 / ul^2 * 2 * sqrt(2.0) * sigma * kv1y
    end
    return real.(FFTW.ifft(dtx)), real.(FFTW.ifft(dty))
end

function _besselk(nu, z)
    try
        return SpecialFunctions.besselk(nu, z)
    catch err
        err isa SpecialFunctions.AmosException || rethrow()
        return zero(z)
    end
end

# Python-style element index: negative positions count from the end.
_py(i::Int, n::Int) = i < 0 ? i + n : i

# Separation bubble on the wind-aligned computational grid (AeoLiS `WindShear.separation`).
# Rows of the Python grid are columns here; indices below are 0-based as in AeoLiS and
# converted with `+ 1` at each access.
function _separation(x::AbstractMatrix, z::AbstractMatrix, dx, sumtaux, c, mu_b,
                     filter_iterations, y_filter)
    nxc, nyc = size(z)
    Z(i, j) = z[_py(i, nxc) + 1, j + 1]
    dzx = zeros(nxc, nyc)
    for j in 0:(nyc - 1)
        for i in 0:(nxc - 3)
            dzx[i + 1, j + 1] = rad2deg(atan((Z(i + 2, j) - Z(i, j)) / (2 * dx)))
        end
        dzx[nxc - 1, j + 1] = dzx[nxc - 2, j + 1]
        dzx[nxc, j + 1] = dzx[nxc - 1, j + 1]
    end
    stall = Float64.((abs.(dzx) .> mu_b) .& (dzx .< 0))
    # numpy.logical_and(a, b, out): the third condition of AeoLiS is its output array
    add = (stall[2:(end - 1), :] .== 0) .& (stall[1:(end - 2), :] .> 0)
    stall[2:(end - 1), :] .+= add
    bubble = zeros(Int, nxc, nyc)
    bubble[1:(end - 1), :] .= (stall[1:(end - 1), :] .== 0) .& (stall[2:end, :] .> 0)
    bubble[2:end, :] .= bubble[1:(end - 1), :]
    bubble[end, :] .= 0
    zsep = zeros(nxc, nyc)
    zsep_new = zeros(nxc, nyc)
    idir = sumtaux >= 0 ? 1 : -1
    for j in 0:(nyc - 1), i in 0:(nxc - 1)    # numpy.where order: by row, then column
        bubble[i + 1, j + 1] == 1 || continue
        start = _py(i + idir * 5, nxc)
        neg = findfirst(>=(0), view(dzx, (start + 1):nxc, j + 1))
        zbrink = neg === nothing ? Z(i, j) : Z(i, j) - Z(i + idir * 5 + idir * (neg - 1), j)
        dzdx0 = (Z(i, j) - Z(i - 3, j)) / (3 * dx)
        a = dzdx0 / c
        ls = min(max(3 * zbrink / (2 * c) * (1 + a / 4 + a^2 / 8), 0.1), 200.0)
        a2 = -3 * zbrink / ls^2 - 2 * dzdx0 / ls
        a3 = 2 * zbrink / ls^3 + dzdx0 / ls^2
        i_max = min(i + trunc(Int, ls / dx) + 1, nxc - 1)
        for k in i:(i_max - 1)
            xs = idir == 1 ? x[k + 1, j + 1] - x[i + 1, j + 1] :
                 -(x[k + 1, j + 1] - x[i + 1, j + 1])
            zsep_new[k + 1, j + 1] = a3 * xs^3 + a2 * xs^2 + dzdx0 * xs + Z(i, j)
        end
        for k in 1:nxc
            zsep[k, j + 1] = max(z[k, j + 1], zsep[k, j + 1], zsep_new[k, j + 1])
        end
        for _ in 1:filter_iterations
            zsep_new = zeros(nxc, nyc)
            Cut = 1.5
            dk = 2.0 * pi / maximum(x)
            kk = 0:(nxc - 1)
            zf = FFTW.fft(complex.(zsep[:, j + 1]))
            zf .*= exp.(-(dk .* kk .* dx) .^ 2 ./ (2 * Cut^2))
            zs = real.(FFTW.ifft(zf))
            Zs(q) = zs[_py(q, nxc) + 1]
            zb1 = neg === nothing ? Zs(i) : Zs(i) - Zs(i + idir * 5 + idir * (neg - 1))
            dzdx1 = (Zs(i) - Zs(i - 3)) / (3 * dx)
            a = dzdx1 / c
            ls = min(max(3 * zb1 / (2 * c) * (1 + a / 4 + a^2 / 8), 0.1), 200.0)
            a2 = -3 * zb1 / ls^2 - 2 * dzdx1 / ls
            a3 = 2 * zb1 / ls^3 + dzdx1 / ls^2
            i_max1 = min(i + idir * trunc(Int, ls / dx), nxc - 1)
            for k in i:(i_max1 - 1)
                xs = idir == 1 ? x[k + 1, j + 1] - x[i + 1, j + 1] :
                     -(x[k + 1, j + 1] - x[i + 1, j + 1])
                zsep_new[k + 1, j + 1] = a3 * xs^3 + a2 * xs^2 + dzdx1 * xs + zb1
            end
            for k in 1:nxc
                zsep[k, j + 1] = max(z[k, j + 1], zsep[k, j + 1], zsep_new[k, j + 1])
            end
        end
    end
    y_filter && (zsep = gaussian_filter1d(zsep, 0.2, 2))
    for k in eachindex(zsep)
        zsep[k] < z[k] && (zsep[k] = z[k])
    end
    return zsep
end
