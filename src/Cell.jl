__precompile__()

module Cell

import SparseArrays: sparse, SparseMatrixCSC
import NamedArrays: NamedArray

CountMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Integer}
NamedCountMatrix{T} = NamedArray{T} where {T <: Integer}

include("utils.jl")
include("input.jl")

export CountMatrix, NamedCountMatrix
export read_10X, read_10X_h5, convert_counts, filter_counts
end # module
