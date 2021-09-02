using Severo
using DataFrames
using HDF5
import CategoricalArrays: CategoricalValue
import Clustering: randindex, counts

include("time_calls.jl")

@time_calls function jl_pipeline(X::NamedCountMatrix, hvf_ref::Vector{String})
    X = filter_counts(X, min_cells=3, min_features=200)
    hvf = find_variable_features(X, 2000, method=:vst)

    Y = normalize_cells(X, scale_factor=1e4, method=:lognormalize, dtype=Float32)

    S = scale_features(Y[:,hvf_ref], scale_max=10)
    em = embedding(S, 150, method=:pca, algorithm=:arpack)
    snn = shared_nearest_neighbours(em, 20, dims=1:50, ntables=50)
    lbls = cluster(snn, resolution=0.5)
    res = umap(em, dims=1:50, metric=:ann)
    @assert eltype(res) == Float32

    de = find_all_markers(Y, lbls; only_pos=true, min_pct=0.25, logfc_threshold=0.25, log=true)

    (names(hvf,1), lbls, res, de)
end

import HDF5: write
function write(parent::Union{HDF5.File, HDF5.Group}, name::AbstractString, df::DataFrame; pv...)
    g = HDF5.create_group(parent, name)

    maybe_convert(x::AbstractVector{T}) where {T} = convert(Vector{T}, x)
    maybe_convert(x::AbstractVector{T}) where {R, T <: CategoricalValue{R}} = convert(Vector{R}, x)

    for n in names(df)
        x = df[!,n]
        x = maybe_convert(x)
        write(g, n, x; pv...)
    end
end

function post_process(f, prefix, x)
  ds = HDF5.create_group(f, prefix)
  t, vals = x

  write(ds, "t/keys", collect(keys(t)))
  write(ds, "t/vals", collect(values(t)))

  hvf = convert(Vector{String}, vals[1])
  write(ds, "hvf", hvf)

  lbls = convert(Vector{Int}, vals[2])
  write(ds, "lbls", lbls)

  umap = convert(Matrix{Float64}, vals[3])
  write(ds, "umap", umap)

  de = convert(DataFrame, vals[4])
  write(ds, "de", de)

  t, hvf, lbls, umap
end

for fname in ARGS
  dataset, _ = splitext(basename(fname))
  try
    h5open("$dataset.h5", "cw") do io
      @assert haskey(io, "R")
      hvf_R = read(io, "R/hvf")

      if !haskey(io, "jl32")
        println("Processing $dataset ($fname)")
        X = read_data(fname)
        t = @elapsed x = jl_pipeline(X, hvf_R)
        if t < 100
          println("short pipeline, running again")
          x = jl_pipeline(X, hvf_R)
        end
        post_process(io, "jl32", x)
      end
    end
  catch e
      println(stderr, "Failed to process $fname")
  end
end
