using Severo
using PyCall
using DataFrames
using HDF5
import CategoricalArrays: CategoricalValue
import Clustering: randindex, counts

include("memusage.jl")
mempeak = memusage()
include("time_calls.jl")

pushfirst!(PyVector(pyimport("sys")."path"), "/shared/common/software/Python/3.8.2-GCCcore-9.3.0/lib/python3.8/site-packages")
sc = pyimport("scanpy")

py"""
def convert_matrix(i, p, x, dim, rn, cn):
  from scipy.sparse import csc_matrix
  from anndata import AnnData
  X = csc_matrix((x, i, p), shape=dim)
  return AnnData(X, {'row_names':rn}, {'col_names':cn}, dtype=X.dtype.name)

def variable_features(X):
   hvf = X.var.highly_variable
   rank = X.var.highly_variable_rank[hvf].sort_values()
   return rank.index

def set_idents(X, lbls):
  import pandas as pd
  X.obs['louvain'] = pd.Categorical(lbls)
  return X

def de_to_df(X):
  import pandas as pd
  de = X.uns['rank_genes_groups']
  frames = [pd.melt(pd.DataFrame.from_records(de[n]), var_name='group', value_name=n).set_index('group', append=True) for n in ["names", "logfoldchanges", "pvals", "pvals_adj"]]
  df = pd.DataFrame.join(frames[0], frames[1:]).reset_index('group')
  df.group = df.group.astype('int')
  return df
"""
function matrix_to_py(X::NamedCountMatrix)
  x = X.array
  py"convert_matrix"(x.rowval .- 1, x.colptr .- 1, convert(Vector{Float32},x.nzval), size(X), names(X,1), names(X, 2))
end

function CreateSeuratObject(X::PyObject; min_cells=3, min_features=200)
  sc.pp.filter_cells(X, min_genes=min_features)
  sc.pp.filter_genes(X, min_cells=min_cells)
  X
end

function NormalizeData(X::PyObject; normalization_method="LogNormalize", scale_factor=1e4)
  @assert normalization_method == "LogNormalize"
  sc.pp.normalize_total(X, target_sum=scale_factor)
  sc.pp.log1p(X)
  X
end

function FindVariableFeatures(X::PyObject; selection_method="vst", nfeatures=2000)
  @assert selection_method == "vst"
  sc.pp.highly_variable_genes(X, flavor="seurat_v3", n_top_genes=nfeatures)
  X
end

VariableFeatures(X::PyObject) = py"variable_features"(X)
function ScaleData(X::PyObject; features=nothing)
  Y = if features !== nothing
    py"$X[:, $features]"
  else
    X
  end

  sc.pp.scale(Y, max_value=10)
  Y
end

function RunPCA(X::PyObject; npcs=50, features=nothing)
  Y = if features !== nothing
    py"$X[:, $features]"
  else
    X
  end

  sc.tl.pca(Y, svd_solver="arpack", n_comps=npcs)
  Y
end

function FindNeighbors(X::PyObject; dims=1:50)
  sc.pp.neighbors(X, n_neighbors=20, n_pcs=maximum(dims))
  X
end

function FindClusters(X::PyObject; resolution=0.8)
  sc.tl.louvain(X, resolution=resolution)
  X
end

function Idents(X::PyObject)
  X.obs.louvain.astype("int")
end

function RunUMAP(X::PyObject; dims=1:50)
  sc.pp.neighbors(X, n_neighbors=30, metric="cosine", n_pcs=maximum(dims))
  sc.tl.umap(X, n_components=2, maxiter=300, min_dist=0.3, init_pos="spectral")
  X
end

function umap_coords(X::PyObject)
  py"$(X).obsm['X_umap']"
end

function pd_to_df(df_pd)
    df= DataFrame()
    for col in df_pd.columns
        df[!, col] = if col == "names"
          convert(Vector{String}, getproperty(df_pd, col).values)
        else
          getproperty(df_pd, col).values
        end
    end
    df
end

function FindAllMarkers(X::PyObject, lbls=nothing; only_pos=true, min_pct=0.25, logfc_threshold=0.25)
  if lbls !== nothing
    py"set_idents"(X, lbls)
  end

  sc.tl.rank_genes_groups(X, "louvain", method="wilcoxon", tie_correct=true)
  df = py"de_to_df"(X)
  pd_to_df(df)
end

@time_calls function py_pipeline(X)
  A = matrix_to_py(X)
  data = CreateSeuratObject(A, min_cells=3, min_features=200)
  data = FindVariableFeatures(data, selection_method="vst", nfeatures=2000)
  hvf = VariableFeatures(data)

  Y = NormalizeData(data, normalization_method="LogNormalize", scale_factor=1e4)

  data = ScaleData(Y, features=hvf)
  data = RunPCA(data, npcs=150, features=hvf)

  data = FindNeighbors(data, dims=1:50)
  data = FindClusters(data, resolution=0.5)
  data = RunUMAP(data, dims=1:50)

  idents = Idents(data)
  de = FindAllMarkers(Y, idents; only_pos=true, min_pct=0.25, logfc_threshold=0.25)

  umap = umap_coords(data)
  hvf, idents, umap, de
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

function post_process(ds, prefix, x)
  t, vals = x

  write(ds, "$(prefix)_t/keys", collect(keys(t)))
  write(ds, "$(prefix)_t/vals", collect(values(t)))

  hvf = convert(Vector{String}, vals[1])
  write(ds, "$(prefix)_hvf", hvf)

  lbls = convert(Vector{Int}, vals[2])
  write(ds, "$(prefix)_lbls", lbls)

  umap = convert(Matrix{Float64}, vals[3])
  write(ds, "$(prefix)_umap", umap)

  de = convert(DataFrame, vals[4])
  write(ds, "$(prefix)_de", de)

  t, hvf, lbls, umap
end

for fname in ARGS
  dataset, _ = splitext(basename(fname))

  try
    h5open("$dataset.h5", "cw") do io
      if !haskey(io, "py_t")
        println("Processing $dataset ($fname)")
        X = read_data(fname)
        post_process(io, "py", py_pipeline(X))
      end
    end
  catch e
      println(stderr, "Failed to process $fname")
  end
end
