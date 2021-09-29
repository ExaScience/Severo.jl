import StatsBase: sample
import Statistics: cor, quantile!
using HDF5
using CSV
using DataFrames
using Distances
import CategoricalArrays: CategoricalArray, cut

impls = ["py", "R", "jl", "jl32"]
nsamples = 10_000
metric = CosineDist()

function utri2vec(U::Matrix{T}) where T
    s = size(U, 1)
    [U[i,j] for i in 1:s for j in i:s]
end

function groupby_boxplot(by::CategoricalArray, val::AbstractVector{S}) where {S, T}
    l = length(by.pool)
    refs = by.refs

    res = zeros(S, (5, l))
    for i in convert(UnitRange{eltype(refs)},1:l)
        sel = refs .== i
        q = view(res, :, i)
        v = view(val, sel)
        quantile!(q, v, [0.0, 0.25, 0.5, 0.75, 1.0])

        iqr = q[4] - q[2]
        if q[1] < q[2] - 1.5 * iqr
            q[1] = (i = findfirst(>(q[2] - 1.5 * iqr), v); v[i])
        end

        if q[5] > q[4] + 1.5 * iqr
            q[5] = (i = findlast(<(q[4] + 1.5 * iqr), v); v[i])
        end
    end

    res
end

df = DataFrame()
for fname in filter(endswith(".h5"), readdir())
    name, _ = splitext(fname)
    println(stderr, "Processing $name")

    it = if startswith(name, "subsample")
        parts = split(name, "_")
        parse(Int64, parts[3])
    else
        1
    end

    h5open(fname, "r") do io
        if haskey(io, "coordinates")
            X = read(io, "coordinates")
            m,n = size(X)
            J = sample(1:m, nsamples)
            X = X[J,:]'

            orig_dists = utri2vec(pairwise(metric, X))
            cuts = cut(orig_dists, 50)

            for impl in impls
                haskey(io, "$impl/umap") || continue

                println(stderr, "Processing $name:$impl")
                umap = read(io, "$impl/umap")[J,:]'
                dists = utri2vec(pairwise(metric, umap))

                r = cor(orig_dists, dists)
                qs = groupby_boxplot(cuts, dists)

                for (i,q) in enumerate(eachcol(qs))
                    q = NamedTuple(zip((:ymin, :lower, :middle, :upper, :ymax), q))
                    push!(df, merge((dataset=name, it=it, size=m, implementation=impl, pearson=r, cut=i), q))
                end
            end
        else
            println(stderr, "Failed to process $name: missing coordinates")
        end
    end
end
CSV.write(stdout, df)
