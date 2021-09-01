using HDF5
using CSV
using DataFrames

const steps_jl = ["filter_counts", "normalize_cells", "find_variable_features", "scale_features", "embedding", "shared_nearest_neighbours", "cluster", "umap", "find_all_markers"]
const steps_R = ["CreateSeuratObject", "NormalizeData", "FindVariableFeatures", "ScaleData", "RunPCA", "FindNeighbors", "FindClusters", "RunUMAP", "FindAllMarkers"]

function read_time(ds, prefix)
    k = read(ds, "$(prefix)_t/keys")
    v = read(ds, "$(prefix)_t/vals")
    replace!(k, "normalize"=>"normalize_cells")
    Dict(zip(k,v))
end

df = DataFrame()
for fname in readdir()
    try
        name, _ = splitext(fname)
        it, size = if name == "fullsample"
            1, 1_000_000
        else
            parts = split(name, "_")
            it = parse(Int64, parts[3])
            size = parse(Int64, parts[5])
            it, size
        end

        h5open(fname, "r") do io
            for (prefix, steps) in [("R", steps_R), ("jl", steps_jl), ("jl32", steps_jl), ("py", steps_R)]
                if haskey(io, "$(prefix)_t")
                    lbls = read(io, "$(prefix)_lbls")
                    t = read_time(io, prefix)
                    vals = [get(t, k, missing) for k in steps]
                    append!(df, DataFrame(dataset=name, it=it, size=size, implementation=prefix, step=steps_R, t=vals, clusters=length(unique(lbls))))
                end
            end
            df
        end
    catch e
        println(stderr, "Failed to process $fname")
        #showerror(stderr, e)
    end
end
CSV.write(stdout, df)
