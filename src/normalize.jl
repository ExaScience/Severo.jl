import SparseArrays: SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange

function row_norm(A::SparseMatrixCSC{T}; scale_factor=1.0) where {T <: Integer}
    B = similar(A, Float64)

    s = sum(A, dims=2)
    @inbounds for i in 1:size(A,2)
        @inbounds for j in nzrange(A, i)
            nonzeros(B)[j] = scale_factor * nonzeros(A)[j] / s[rowvals(A)[j]]
        end
    end

    B
end

function log_norm(A::SparseMatrixCSC{T}; scale_factor=1.0) where {T <: Integer}
    B = row_norm(A, scale_factor=scale_factor)
    nonzeros(B) .= log1p.(nonzeros(B))
    B
end

"""
    normalize(X::NamedCountMatrix; method=:lognormalize, scale_factor=1.0)

Normalize count data with different methods:

    - `lognormalize`: feature counts are divided by the total count per cell, scaled by `scale_factor` and then log1p transformed.
    - `relativecounts`: feature counts are divided by the total count per cell and scaled by `scale_factor`.

**Arguments**:

    - `X`: the labelled count matrix to normalize
    - `method`: normalization method to apply
    - `scale_factor`: the scaling factor

**Return values**:

A labelled data matrix
"""
function normalize(X::NamedCountMatrix; method=:lognormalize, scale_factor=1e4)
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    if method == :lognormalize
        f = log_norm
    elseif method == :relativecounts
        f = row_norm
    else
        error("unknown normalization method: $method")
    end

    S = f(X.array; scale_factor=scale_factor)
    NamedArray(S, X.dicts, X.dimnames)
end
