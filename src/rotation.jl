# Rotated computational grids (AeoLiS `rotation.rotationClass`) for processes computed along
# wind-aligned transects: Okin's vegetation shear, and the `quasi2d` (Kroy et al. 2002) and
# `duna2d` shear perturbations with their separation bubble.

"""
    RotationGrid(dx, dy, buffer_width)

A wind-aligned computational grid with cells `dx` by `dy` and a buffer around the input grid.
"""
mutable struct RotationGrid
    dx::Float64
    dy::Float64
    buffer_width::Float64
    x0::Float64
    y0::Float64
end
RotationGrid(dx, dy, bw) = RotationGrid(dx, dy, bw, 0.0, 0.0)

# The grid as in `WindShear`; the interpolation is the same as well.
function _rotation_setup(rg::RotationGrid, s::AeolisState)
    udir = s.udir[1, 1]
    u_angle = 270.0 - udir
    xi, yi = _computational_grid!(rg, s.x, s.y, udir)
    o = (rg.x0, rg.y0)
    gcx, gcy = rotate(xi, yi, -u_angle; origin=o)
    gcz = _shear_interpolate(rg, s.x, s.y, s.zb, gcx, gcy, 0.0)
    return u_angle, gcx, gcy, gcz, o
end

"""
    okin_shear!(rg, s, m, c1, intercept)

The shear velocity reduced in the lee of vegetation (Okin 2008) along wind-aligned rows,
interpolated back to the model grid.
"""
function okin_shear!(rg::RotationGrid, s::AeolisState, m::AeolisModel, c1, intercept)
    u_angle, gcx, gcy, _, o = _rotation_setup(rg, s)
    gcx, gcy = rotate(gcx, gcy, u_angle; origin=o)
    gix, giy = rotate(s.x, s.y, u_angle; origin=o)
    gus = _shear_interpolate(rg, gix, giy, s.ustar, gcx, gcy, 0.0)
    ghv = _shear_interpolate(rg, gix, giy, s.hveg, gcx, gcy, 0.0)
    gus = _okin_rows(gcx, gus, ghv, c1, intercept)
    return _shear_interpolate(rg, gcx, gcy, gus, gix, giy, 0.0)
end

"""
    rotated_shear!(rg, s, m)

Shear stress perturbed along wind-aligned transects by `method_shear` (`:duna2d`, the
Duna integral; `:quasi2d`, Kroy et al. 2002), with the quasi-2-D separation bubble. Returns
`(tau, zsep, hsep)` on the model grid.
"""
function rotated_shear!(rg::RotationGrid, s::AeolisState, m::AeolisModel)
    p = m.p
    u_angle, gcx, gcy, gcz, o = _rotation_setup(rg, s)
    gcx, gcy = rotate(gcx, gcy, u_angle; origin=o)
    gix, giy = rotate(s.x, s.y, u_angle; origin=o)
    gtau = _shear_interpolate(rg, gix, giy, s.tau, gcx, gcy, 0.0)
    if p.method_shear == :duna2d
        gtau = _shear_duna(gcx, gcz, gtau, 3.0, 1.0)
    else
        z0 = roughness_length(s, m)
        gtau = _shear_kroy(gcx, gcz, gtau, p.L, p.l, p.kappa, p.rhoa * 1.5e-5, p.rhoa, z0)
    end
    zsep, hsep, gtau = _separation_quasi2d(gcx, gcz, gtau, rg.dx, p)
    zs = _shear_interpolate(rg, gcx, gcy, zsep, gix, giy, 0.0)
    hs = _shear_interpolate(rg, gcx, gcy, hsep, gix, giy, 0.0)
    tau = _shear_interpolate(rg, gcx, gcy, gtau, gix, giy, 0.0)
    return tau, zs, hs
end

# AeoLiS `compute_shear_perturbation_1D`: the Duna integral along each row.
function _shear_duna(x, z, tau, alfa, beta)
    nxc, nyc = size(x)
    dx = abs(x[2, 1] - x[1, 1])
    out = copy(tau)
    n = nxc - 1
    for j in 1:nyc
        dz = zeros(nxc)
        dz[2:(end - 1)] .= (z[3:end, j] .- z[1:(end - 2), j]) ./ 2 ./ dx
        for i in 0:n
            integ = 0.0
            for jj in (i - n):(i - 1)
                jj != 0 && (integ += dz[i - jj + 1] / (jj * pi))
            end
            r = max(alfa * (integ + beta * dz[i + 1]) + 1, 0.1)
            out[i + 1, j] = tau[i + 1, j] * r
        end
    end
    return out
end

# AeoLiS `compute_shear_perturbation_kroy1D` with `Kroy2002.calc_tau_p`.
function _shear_kroy(x, z, tau, dune_length, l, kappa, mu, rho, z0)
    nxc, nyc = size(x)
    xv = x[:, 1]
    dx = abs(xv[2] - xv[1])
    x2 = xv[2:end] .- dx / 2
    out = copy(tau)
    out[isnan.(out)] .= 0
    for j in 1:nyc
        tp = _kroy_tau_p(xv, z[:, j], z0, dune_length, kappa)
        tpi = length(tp) == length(x2) ? [clamped_interp(xi, x2, tp) for xi in xv] : tp
        tpi[isnan.(tpi)] .= 0
        tk = tau[:, j] .* (1 .+ tpi)
        for k in eachindex(tk)
            isinf(tk[k]) && (tk[k] = tau[k, j])
        end
        out[:, j] .= tk
    end
    for k in eachindex(out)
        (isinf(out[k]) || isnan(out[k]) || out[k] < 0) && (out[k] = tau[k])
    end
    return out
end

function _kroy_tau_p(x, z, z0, dune_length, kappa)
    L = 0.5 * dune_length
    Phi = L / z0
    f(v) = v - 2.0 * kappa^2 * Phi / log(v)
    df(v) = 1 + 2.0 * kappa^2 * Phi / (v * log(v)^2)
    phi = 2.0 * kappa^2 * Phi / log(2.0 * kappa^2 * Phi)
    for _ in 1:50
        step = f(phi) / df(phi)
        phi -= step
        abs(step) < 1.48e-8 && break
    end
    fact = 1 + log(phi) + 2.0 * log(0.5 * pi) + 4.0 * Base.MathConstants.eulergamma
    a = ((log(Phi^2 / log(Phi))^2) / (2.0 * log(phi)^3)) * fact
    b = pi / fact
    n = length(z)
    F = FFTW.rfft(z)
    freq = 2.0 * pi .* (0:(length(F) - 1)) ./ (n * (x[2] - x[1]))
    Ft = a .* (freq .+ 1im .* b .* abs.(freq)) .* F
    return FFTW.irfft(Ft, 2 * (length(F) - 1))
end

# AeoLiS `separation_quasi2d`: the separation bubble on the wind-aligned grid and the reduced
# stress under it.
function _separation_quasi2d(x, z, tau, dx, p)
    if !p.process_separation
        return zeros(size(z)), zeros(size(z)), tau
    end
    zsep = _separation(x, z, dx, 1.0, p.c_b, p.mu_b, p.sep_filter_iterations,
                       p.zsep_y_filter)
    hsep = zsep .- z
    delta = 1.0 / (0.2 * 0.5)
    f = min.(max.(1.0 .- delta .* hsep, 0.0), 1.0)
    return zsep, hsep, tau .* f
end
