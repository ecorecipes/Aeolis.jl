# Sand fences (AeoLiS `fences.py`): fences are buried as the bed rises, and reduce the shear
# velocity in their lee (Okin 2008) on a wind-aligned computational grid.

"""
    initialize_fences!(s, m)

Fence heights from `fence_file` and the level of their tops above the initial bed.
"""
function initialize_fences!(s::AeolisState, m::AeolisModel)
    p = m.p
    p.process_fences || return s
    s.fence_height .= _as_grid(eltype(s.zb), p.fence_file, m.nx + 1, m.ny + 1, :fence_file)
    s.fence_base .= s.zb
    s.fence_top .= s.fence_base .+ s.fence_height
    s.fence_height_init .= s.fence_height
    return s
end

"""
    update_fences!(s, m)

The exposed fence heights (fences lower than 0.1 m or higher than 1.5 m are ignored), and the
shear velocity reduced in their lee.
"""
function update_fences!(s::AeolisState, m::AeolisModel)
    _update_fence_height!(s)
    _fence_shear!(s, m)
    velocity_stress!(s, m)
    return s
end

function _update_fence_height!(s)
    s.fence_height .= s.fence_top .- s.zb
    for k in eachindex(s.fence_height)
        s.fence_height_init[k] < 0.1 && (s.fence_height[k] = 0)
        s.fence_top[k] < 0.1 && (s.fence_height[k] = 0)
        s.fence_height[k] < 0.1 && (s.fence_height[k] = 0)
        if s.fence_height[k] > 1.5
            s.fence_height[k] = 0
            s.fence_height_init[k] = 0
        end
    end
    return s
end

# RegularGridInterpolator over the axis-aligned grid (x, y), fill 0 (AeoLiS
# `fences.interpolate`).
_fence_interp(x, y, z, xi, yi) = regular_grid_interp(x[:, 1], y[1, :], z, xi, yi, 0.0)

# AeoLiS `fence_shear2d` (used for transects as well, which are quasi-2-D).
function _fence_shear!(s::AeolisState, m::AeolisModel)
    p = m.p
    x, y, zf = s.x, s.y, s.fence_height
    dx = max(p.dx / 2, 0.25)
    dy = dx
    udir = s.udir[1, 1]
    udir < 0 && (udir += 360)
    udir = udir - floor(udir / 360) * 360
    udir in (0.0, 360.0, -360.0, -180.0, 180.0) && (udir += 0.00001)
    x0, y0 = mean(x), mean(y)
    D = sqrt((maximum(x) - minimum(x))^2 + (maximum(y) - minimum(y))^2) + 2 * 100.0
    xi, yi = exact_grid(x0 - D / 2, x0 + D / 2, y0 - D / 2, y0 + D / 2, dx, dy)
    alpha = udir + 90.0
    o = (x0, y0)
    xc, yc = rotate(xi, yi, alpha; origin=o)
    zfc = _fence_interp(x, y, zf, xc, yc)
    usx = _fence_interp(x, y, s.ustars, xc, yc)
    # Okin reduction along the rows of the computational grid
    nxc, nyc = size(zfc)
    mult_all = ones(nxc, nyc)
    c1, intercept = p.okin_c1_fence, p.okin_initialred_fence
    for j in 1:nyc
        red_all = zeros(nxc)
        for g in 1:nxc
            h = zfc[g, j]
            h > 0 || continue
            for k in 1:nxc
                xrel = (udir >= 180 && udir <= 360) ? xc[k, j] - xc[g, j] :
                       -(xc[k, j] - xc[g, j])
                mult = (xrel >= 0 && xrel / h < 20) ?
                       intercept + (1 - intercept) * (1 - exp(-xrel * c1 / h)) : 1.0
                red = 1 - mult
                red < 0.00001 && (red = 0.0)
                red > 1 && (red = 1.0)
                xrel < 0 && (red = 0.0)
                red_all[k] += red
            end
        end
        mult_all[:, j] .= 1 .- min.(red_all, 1.0)
    end
    mult_all[1, :] .= 1
    mult_all[end, :] .= 1
    mult_all[:, 1] .= 1
    mult_all[:, end] .= 1
    gux = usx .* mult_all
    guy = zeros(size(gux))
    gux, guy = rotate(gux, guy, alpha)
    gix, giy = rotate(x, y, -alpha; origin=o)
    xc, yc = rotate(xc, yc, -alpha; origin=o)
    guy, gux = rotate(gux, guy, -alpha)
    ux = _fence_interp(xc, yc, gux, gix, giy)
    uy = _fence_interp(xc, yc, guy, gix, giy)
    ux, uy = rotate(ux, uy, alpha)
    for a in (ux, uy)
        a[:, 1] .= a[:, 2]
        a[1, :] .= a[2, :]
        a[:, end] .= a[:, end - 1]
        a[end, :] .= a[end - 1, :]
    end
    s.ustars .= ux
    s.ustarn .= uy
    s.ustar .= sqrt.(s.ustars .^ 2 .+ s.ustarn .^ 2)
    return s
end
