# copyright imec - evaluation license - not for distribution

__precompile__()

module Severo

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 3
end

import SparseArrays: sparse, SparseMatrixCSC, SparseVector, SparseColumnView
import NamedArrays: NamedArray, NamedMatrix, NamedVector, names, dimnames

const CountMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Integer}
const NamedCountMatrix{T} = NamedArray{T} where {T <: Integer}

const DataMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Real}
const NamedDataMatrix{T} = NamedArray{T} where {T <: Real}

const NeighbourGraph{T} = NamedArray{T, 2, SparseMatrixCSC{T, Int64}} where {T <: Real}

const SparseVec = Union{SparseColumnView, SparseVector}

import Requires: @require

try
    include(joinpath(dirname(@__DIR__), "deps","deps.jl"))
catch e
    error("Severo.jl not properly configured, please run `Pkg.build(\"Severo\")`.")
end

include("utils.jl")
include("input.jl")
include("filtering.jl")
include("normalize.jl")
include("scaling.jl")
include("variablefeatures.jl")
include("embedding.jl")
include("neighbours.jl")
include("clustering.jl")
include("datasets.jl")
include("diffexpr.jl")
include("printing.jl")

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("visualization.jl")
end

export CountMatrix, NamedCountMatrix, DataMatrix, NamedDataMatrix
export read_10X, read_10X_h5, read_h5, read_h5ad, write_h5ad, read_csv, read_geo, read_loom, read_data
export convert_counts, filter_counts, filter_features, filter_cells, dataset
export normalize_cells, scale_features, find_variable_features
export embedding, pca, umap
export nearest_neighbours, shared_nearest_neighbours, jaccard_index
export cluster, prefilter_markers, find_markers, filter_rank_markers, find_all_markers
end # module