module AeolisNCDatasetsExt

using Aeolis
using NCDatasets
using Dates

# netCDF dimension names of AeoLiS output for the dimensions of the state table
const NCDIMS = Dict(:nx => "s", :ny => "n", :nlayers => "layers",
                    :nfractions => "fractions")

function Aeolis.write_netcdf(file::AbstractString, sol::Aeolis.AeolisSolution)
    s = sol.u
    p = sol.prob.p
    nxp, nyp = size(s.zb)
    transect = nyp == 3
    nn = transect ? 1 : nyp
    rows = transect ? (2:2) : (1:nyp)
    nt = length(sol.t)
    NCDataset(file, "c") do ds
        ds.attrib["Conventions"] = "CF-1.6"
        ds.attrib["source"] = "Aeolis.jl"
        ds.attrib["date_created"] = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MMZ")
        defDim(ds, "s", nxp)
        defDim(ds, "n", nn)
        defDim(ds, "time", Inf)
        defDim(ds, "layers", size(s.mass, 3))
        defDim(ds, "fractions", size(s.mass, 4))
        defVar(ds, "s", Float32.(0:(nxp - 1)), ("s",);
               attrib=["long_name" => "s-coordinate"])
        defVar(ds, "n", Float32.(0:(nn - 1)), ("n",);
               attrib=["long_name" => "n-coordinate"])
        defVar(ds, "x", Float32.(s.x[:, rows]), ("s", "n");
               attrib=["long_name" => "x-coordinate", "units" => "m"])
        defVar(ds, "y", Float32.(s.y[:, rows]), ("s", "n");
               attrib=["long_name" => "y-coordinate", "units" => "m"])
        defVar(ds, "layers", Float32.(0:(size(s.mass, 3) - 1)), ("layers",))
        defVar(ds, "fractions", Float32.(p.grain_size), ("fractions",);
               attrib=["long_name" => "sediment fractions", "units" => "m"])
        tv = defVar(ds, "time", Float64, ("time",);
                    attrib=["units" => "seconds since $(p.refdate)",
                            "calendar" => "julian"])
        tv[1:nt] = sol.t
        for (key, a) in sol.saved
            base = Symbol(first(split(String(key), r"_(avg|sum|var|min|max)$")))
            info = Aeolis.state_info(base)
            dims = info.dims
            sel = selectdim(a, 2, rows)
            # (x, y, extra..., time) -> (extra reversed..., x, y, time) for C-order (time, n, s, ...)
            nd = ndims(a)
            perm = vcat(collect((nd - 1):-1:3), [1, 2, nd])
            data = Float32.(permutedims(sel, perm))
            ncd = vcat([NCDIMS[d] for d in reverse(dims[3:end])], ["s", "n", "time"])
            defVar(ds, String(key), data, Tuple(ncd);
                   attrib=["long_name" => String(key),
                           "units" => info.unit == "-" ? "1" : info.unit])
        end
    end
    return file
end

end
