using HDF5
using CSV
using DataFrames

const steps_jl = ["filter_counts", "normalize", "find_variable_features", "scale_features", "embedding", "shared_nearest_neighbours", "cluster", "umap", "find_all_markers"]
const steps_R = ["CreateSeuratObject", "NormalizeData", "FindVariableFeatures", "ScaleData", "RunPCA", "FindNeighbors", "FindClusters", "RunUMAP", "FindAllMarkers"]

function read_time(ds, prefix)
    k = read(ds, "$(prefix)_t/keys")
    v = read(ds, "$(prefix)_t/vals")
    Dict(zip(k,v))
end

df = DataFrame()
for fname in readdir()
    try
        name, _ = splitext(fname)
        parts = split(name, "_")
        it = parse(Int64, parts[3])
        size = parse(Int64, parts[5])

        h5open(fname, "r") do io
            for (prefix, steps) in [("R", steps_R), ("jl", steps_jl), ("py", steps_R)]
                if haskey(io, "$(prefix)_t")
                    t = read_time(io, prefix)
                    vals = [getindex(t, k) for k in steps]
                    append!(df, DataFrame(dataset=name, it=it, size=size, implementation=prefix, step=steps_R, t=vals))
                end
            end
            df
        end
    catch
        println("Failed to process $fname")
    end
end
CSV.write(stdout, df)
