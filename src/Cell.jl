module Cell

import SparseArrays: sparse, SparseMatrixCSC
import NamedArrays: NamedArray

CountMatrix{T} = SparseMatrixCSC{T, Int64} where {T <: Signed}
NamedCountMatrix{T} = NamedArray{T} where {T <: Signed}

include("utils.jl")
include("input.jl")

export CountMatrix, NamedCountMatrix
export read_10X, read_10X_h5, convert_counts, filter_counts
end # module
