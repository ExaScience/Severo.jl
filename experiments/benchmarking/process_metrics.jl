using HDF5
using CSV
using DataFrames
import Clustering: randindex

jaccard(A, B) = (A = Set(A); B = Set(B); length(intersect(A,B)) / length(union(A,B)))

function read_peakmem(ds)
    if haskey(ds, "peakmem")
        read(ds, "peakmem")/1024^2
    else
        missing
    end
end

df = DataFrame(dataset=String[], it=Int64[], size=Int64[], implementation=String[],
     ari=Float64[], ri=Float64[], jaccard=Float64[], peakmem=Union{Missing, Float64}[])
for fname in filter(endswith("h5"), readdir())
    try
        name, _ = splitext(fname)
        it = if startswith(name, "subsample")
            parts = split(name, "_")
            parse(Int64, parts[3])
        else
            1
        end

        h5open(fname, "r") do io
            if haskey(io, "R")
                lbls_ref = read(io, "R/lbls")
                hvf_ref = read(io, "R/hvf")
                size = length(lbls_ref)
                peakmem = read_peakmem(io["R"])

                push!(df, (dataset=name, it=it, size=size, implementation="R", ari=1., ri=1., jaccard=1., peakmem=peakmem))

                for prefix in filter(x->haskey(io, x), ["jl", "jl32", "py"])
                    lbls = read(io, "$(prefix)/lbls")
                    hvf = read(io, "$(prefix)/hvf")
                    ari, ri, _ = randindex(lbls_ref, lbls)
                    j = jaccard(hvf_ref, hvf)
                    peakmem = read_peakmem(io[prefix])

                    push!(df, (dataset=name, it=it, size=size, implementation=prefix, ari=ari, ri=ri, jaccard=j, peakmem=peakmem))
                end

            end
            df
        end
    catch e
        println(stderr, "Failed to process $fname: $e")
    end
end
CSV.write(stdout, df)
