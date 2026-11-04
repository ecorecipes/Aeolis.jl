# Plot recipes (RecipesBase): profiles of a saved variable along a transect over time, or a
# map of a two-dimensional grid at one output time.

RecipesBase.@recipe function f(sol::AeolisSolution, name::Symbol; times=nothing,
                               index=nothing,
                               fraction=1)
    a = sol[name]
    nt = length(sol.t)
    x = sol.u.x
    if size(a, 2) == 3
        # transect: the middle row, at several times
        ks = times === nothing ? unique(round.(Int, range(1, nt; length=min(nt, 5)))) :
             [argmin(abs.(sol.t .- t)) for t in times]
        xlabel --> "x (m)"
        ylabel --> string(name, " (", state_info(name).unit, ")")
        for k in ks
            v = ndims(a) == 3 ? a[:, 2, k] : a[:, 2, fraction, k]
            RecipesBase.@series begin
                label --> string("t = ", round(sol.t[k] / 3600; digits=1), " h")
                x[:, 2], v
            end
        end
    else
        k = index === nothing ? nt : index
        v = ndims(a) == 3 ? a[:, :, k] : a[:, :, fraction, k]
        seriestype --> :heatmap
        aspect_ratio --> :equal
        xlabel --> "x (m)"
        ylabel --> "y (m)"
        title --> string(name, " at t = ", round(sol.t[k] / 3600; digits=1), " h")
        x[:, 1], sol.u.y[1, :], permutedims(v)
    end
end
