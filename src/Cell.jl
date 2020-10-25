__precompile__()

module Cell

import SparseArrays: sparse, SparseMatrixCSC
import NamedArrays: NamedArray

CountMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Integer}
NamedCountMatrix{T} = NamedArray{T} where {T <: Integer}

DataMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Integer}
NamedDataMatrix{T} = NamedArray{T} where {T <: Integer}

try
    include(joinpath(dirname(@__DIR__), "deps","deps.jl"))
catch e
    error("Cell.jl not properly configured, please run `Pkg.build(\"Cell\")`.")
end

include("utils.jl")
include("input.jl")
include("normalize.jl")

export CountMatrix, NamedCountMatrix, DataMatrix, NamedDataMatrix
export read_10X, read_10X_h5, convert_counts, filter_counts
export normalize
end # module
