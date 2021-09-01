using Severo
using RCall
using DataFrames
using HDF5
import CategoricalArrays: CategoricalValue
import Clustering: randindex, counts

include("time_calls.jl")
@rimport Seurat

function matrix_to_R(X::NamedCountMatrix)
  x = X.array
	R"""
(function(i, p, x, dim, rn, cn) {
  A <- new("dgCMatrix")
  A@i <- i
  A@p <- p
  A@x <- x
  A@Dim <- dim
  rownames(A) <- rn
  colnames(A) <- cn
  Matrix::t(A)
})($(x.rowval .- 1), $(x.colptr .- 1), $(convert(Vector{Float64},x.nzval)), $(collect(size(x))), $(names(X,1)), $(names(X, 2)))
	 """
end

umap_coords = R"""
function(data) data[["umap"]]@cell.embeddings
"""

CreateSeuratObject(A::RObject; min_cells=3, min_features=200) = R"Seurat::CreateSeuratObject(counts=$A, min.cells=$min_cells, min.features=$min_features)"
NormalizeData(data::RObject; normalization_method="LogNormalize", scale_factor=1e4) = R"Seurat::NormalizeData(object=$data, normalization.method=$normalization_method, scale.factor=$scale_factor)"
FindVariableFeatures(data::RObject; selection_method="vst", nfeatures=2000) = R"Seurat::FindVariableFeatures(object=$data, selection.method=$selection_method, nfeatures=$nfeatures)"
VariableFeatures(data::RObject) = R"Seurat::VariableFeatures(object=$data)"
Idents(data::RObject) = R"Seurat::Idents(object=$data)"
ScaleData(data::RObject; features=nothing) = R"Seurat::ScaleData(object=$data, features=$features)"
RunPCA(data::RObject; npcs=50, features=nothing) = R"Seurat::RunPCA(object=$data, npcs=$npcs, features=$features)"
FindNeighbors(data::RObject; dims=1:50) = R"Seurat::FindNeighbors(object=$data, dims=$dims)"
FindClusters(data::RObject; resolution=0.8) = R"Seurat::FindClusters(object=$data, resolution=$resolution)"
RunUMAP(data::RObject; dims=1:50) = R"Seurat::RunUMAP(object=$data, dims=$dims)"
FindAllMarkers(data::RObject; only_pos=true, min_pct=0.25, logfc_threshold=0.25) = R"Seurat::FindAllMarkers($data, only.pos=$only_pos, min.pct=$min_pct, logfc.threshold=$logfc_threshold)"

@time_calls function R_pipeline(A)
	data = CreateSeuratObject(A, min_cells=3, min_features=200)
  println(data)
	data = NormalizeData(data, normalization_method="LogNormalize", scale_factor=1e4)
	data = FindVariableFeatures(data, selection_method="vst", nfeatures=2000)

	hvf = VariableFeatures(data)
	data = ScaleData(data, features=hvf)
	data = RunPCA(data, npcs=150, features=hvf)

	data = FindNeighbors(data, dims=1:50)
	data = FindClusters(data, resolution=0.5)
	data = RunUMAP(data, dims=1:50)

  de = FindAllMarkers(data; only_pos=true, min_pct=0.25, logfc_threshold=0.25)

	idents = Idents(data)
	umap = umap_coords(data)
	hvf, idents, umap, de
end

@time_calls function jl_pipeline(X::NamedCountMatrix, hvf_ref::Vector{String})
	X = filter_counts(X, min_cells=3, min_features=200)
	hvf = find_variable_features(X, 2000, method=:vst)

	Y = normalize_cells(X, scale_factor=1e4, method=:lognormalize)

	S = scale_features(Y[:,hvf_ref], scale_max=10)
	em = embedding(S, 150, method=:pca, algorithm=:irlba)
	snn = shared_nearest_neighbours(em, 20, dims=1:50, ntables=50)
	lbls = cluster(snn, resolution=0.5)
	res = umap(em, dims=1:50, metric=:ann)

  de = find_all_markers(Y, lbls; only_pos=true, min_pct=0.25, logfc_threshold=0.25, log=true)

	(names(hvf,1), lbls, res, de, snn, view(em.coordinates, :, 1:50))
end

function modularity(snn, lbls, resolution=0.5)
  n = Severo.Network(snn)
  c = Severo.Clustering(n, lbls, resolution)
  Severo.modularity(c)
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

function read_process(f, prefix)
  println("reading $prefix")
  ds = f[prefix]

  t = begin
    k = read(ds, "t/keys")
    v = read(ds, "t/vals")
    Dict(zip(k,v))
  end

  hvf = read(ds, "hvf")
  lbls = read(ds, "lbls")
  umap = read(ds, "umap")
  t, hvf, lbls, umap
end

for fname in ARGS
  dataset, _ = splitext(basename(fname))
  println("Processing $dataset ($fname)")

  X = read_data(fname)
  h5open("$dataset.h5", "cw") do io
    t_R, hvf_R, lbls_R, umap_R = if !haskey(io, "R")
      A = matrix_to_R(X)
      rv = post_process(io, "R", R_pipeline(A))
      A = nothing;
      GC.gc(); R"gc(reset=T, full=T)"
      rv
    else
      read_process(io, "R")
    end

    if !haskey(io, "jl")
      t = @elapsed x = jl_pipeline(X, hvf_R)
      if t < 100
        println("short pipeline, running again")
        x = jl_pipeline(X, hvf_R)
      end
      t_jl, hvf_jl, lbls_jl, umap_jl = post_process(io, "jl", x)

      snn = last(x)[5].array
      modularity_R = modularity(snn, lbls_R, 0.5)
      modularity_jl = modularity(snn, lbls_jl, 0.5)

      ari, _ = randindex(lbls_R, lbls_jl)
      jaccard(A, B) = (A = Set(A); B = Set(B); length(intersect(A,B)) / length(union(A,B)))
      j = jaccard(hvf_R, hvf_jl)

      println("$dataset, $ari, $j, $modularity_R, $modularity_jl")
      write(io, "comparison", [ari, j, modularity_R, modularity_jl])
      write(io, "coordinates", copy(last(x)[6]))
    end
  end
end
