import SparseArrays: SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange

function log_norm(A::SparseMatrixCSC{T}; scale_factor=1.0) where {T <: Integer}
    B = similar(A, Float64)

    s = sum(A, dims=2)
    for (a, b) in zip(eachcol(A), eachcol(B))
        @inbounds for (i, idx) in enumerate(nonzeroinds(a))
            nonzeros(b)[i] = log1p(scale_factor * nonzeros(a)[i] / s[idx])
        end
#    nonzeros(b) .= log1p.(scale_factor * nonzeros(a) ./ s[nonzeroinds(a)])
    end

    B
end

function row_norm(A::SparseMatrixCSC{T}; scale_factor=1.0) where {T <: Integer}
    B = similar(A, Float64)

    s = sum(A, dims=2)
    for (a, b) in zip(eachcol(A), eachcol(B))
        @inbounds for (i, idx) in enumerate(nonzeroinds(a))
            nonzeros(b)[i] = scale_factor * nonzeros(a)[i] / s[idx]
        end
    end

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
function normalize(X::NamedCountMatrix; method=:lognormalize, scale_factor=1.0)
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
