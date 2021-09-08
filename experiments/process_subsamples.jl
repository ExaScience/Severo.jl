using HDF5
using CSV
using DataFrames

const steps_jl = ["filter_counts", "normalize_cells", "find_variable_features", "scale_features", "embedding", "shared_nearest_neighbours", "cluster", "umap", "find_all_markers"]
const steps_R = ["CreateSeuratObject", "NormalizeData", "FindVariableFeatures", "ScaleData", "RunPCA", "FindNeighbors", "FindClusters", "RunUMAP", "FindAllMarkers"]

function read_time(ds)
    k = read(ds, "t/keys")
    v = read(ds, "t/vals")
    replace!(k, "normalize"=>"normalize_cells")
    Dict(zip(k,v))
end

function find_cellcount(f)
    for ds in f
        if haskey(ds, "lbls")
            return length(ds["lbls"])
        end
    end

    error("cell count not found")
end

df = DataFrame()
for fname in filter(endswith(".h5"), readdir())
    try
        name, _ = splitext(fname)
        h5open(fname, "r") do io
            it = if startswith(name, "subsample")
                parts = split(name, "_")
                it = parse(Int64, parts[3])
                size = parse(Int64, parts[5])
                it
            else
                1
            end

            size = find_cellcount(io)

            for (prefix, steps) in [("R", steps_R), ("jl", steps_jl), ("jl32", steps_jl), ("py", steps_R), ("jl32_opt", steps_jl), ("jl_opt", steps_jl)]
                if haskey(io, prefix)
                    ds = io[prefix]
                    lbls = read(ds, "lbls")
                    t = read_time(ds)
                    vals = [get(t, k, missing) for k in steps]
                    append!(df, DataFrame(dataset=name, it=it, size=size, implementation=prefix, step=steps_R, time=vals, clusters=length(unique(lbls))))
                end
            end
            df
        end
    catch e
        println(stderr, "Failed to process $fname")
        showerror(stderr, e)
    end
end
CSV.write(stdout, df)
