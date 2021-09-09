using HDF5
using CSV
using DataFrames
import Clustering: randindex

jaccard(A, B) = (A = Set(A); B = Set(B); length(intersect(A,B)) / length(union(A,B)))

df = DataFrame()
for fname in filter(endswith("h5"), readdir())
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
            if haskey(io, "R")
                lbls_ref = read(io, "R/lbls")
                hvf_ref = read(io, "R/hvf")

                for prefix in filter(x->haskey(io, x), ["jl", "jl32", "py"])
                    lbls = read(io, "$(prefix)/lbls")
                    hvf = read(io, "$(prefix)/hvf")
                    ari, ri, _ = randindex(lbls_ref, lbls)
                    j = jaccard(hvf_ref, hvf)
                    push!(df, (dataset=name, it=it, size=size, implementation=prefix, ari=ari, ri=ri, jaccard=j))
                end

            end
            df
        end
    catch e
        println(stderr, "Failed to process $fname")
    end
end
CSV.write(stdout, df)
