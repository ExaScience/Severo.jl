__precompile__()

module Cell

import SparseArrays: sparse, SparseMatrixCSC
import NamedArrays: NamedArray, NamedMatrix, NamedVector, names, dimnames

CountMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Integer}
NamedCountMatrix{T} = NamedArray{T} where {T <: Integer}

DataMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Real}
NamedDataMatrix{T} = NamedArray{T} where {T <: Real}

NeighbourGraph{T} = NamedArray{T, 2, SparseMatrixCSC{T, Int64}} where {T <: Real}

try
    include(joinpath(dirname(@__DIR__), "deps","deps.jl"))
catch e
    error("Cell.jl not properly configured, please run `Pkg.build(\"Cell\")`.")
end

include("utils.jl")
include("input.jl")
include("normalize.jl")
include("scaling.jl")
include("variablefeatures.jl")
include("embedding.jl")
include("neighbours.jl")
include("clustering.jl")

export CountMatrix, NamedCountMatrix, DataMatrix, NamedDataMatrix
export read_10X, read_10X_h5, read_h5, read_h5ad
export convert_counts, filter_counts
export normalize, scale, find_variable_features
export embedding, pca, umap
export nearest_neigbours, shared_nearest_neigbours
export cluster
end # module
