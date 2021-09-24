# Severo: a software package for analysis and exploration of single-cell RNA-seq datasets.
# Copyright (c) 2021 imec vzw.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, and Additional Terms
# (see below).

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.

import SparseArrays: SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange

function row_norm(A::SparseMatrixCSC{T}, scale_factor::R=1.0) where {T <: Integer, R <: AbstractFloat}
    B = similar(A, R)

    s = sum(A, dims=2)
    @inbounds for i in 1:size(A,2)
        @inbounds for j in nzrange(A, i)
            nonzeros(B)[j] = scale_factor * nonzeros(A)[j] / s[rowvals(A)[j]]
        end
    end

    B
end

function log_norm(A::SparseMatrixCSC{T}, scale_factor::R=1.0) where {T <: Integer, R <: AbstractFloat}
    B = row_norm(A, scale_factor)
    nonzeros(B) .= log1p.(nonzeros(B))
    B
end

"""
    normalize_cells(X::NamedCountMatrix; method=:lognormalize, scale_factor=1.0)

Normalize count data with different methods:

    - `lognormalize`: feature counts are divided by the total count per cell, scaled by `scale_factor` and then log1p transformed.
    - `relativecounts`: feature counts are divided by the total count per cell and scaled by `scale_factor`.

**Arguments**:

    - `X`: the labelled count matrix to normalize
    - `method`: normalization method to apply
    - `scale_factor`: the scaling factor
    - `dtype`: datatype to be used for the output

**Return values**:

A labelled data matrix
"""
@partial function normalize_cells(X::NamedCountMatrix; method=:lognormalize, scale_factor::Real=1., dtype::Type{T}=Float64) where {T <: AbstractFloat}
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    f = if method == :lognormalize
        log_norm
    elseif method == :relativecounts
        row_norm
    else
        error("unknown normalization method: $method")
    end

    scale_factor = convert(dtype, scale_factor)
    S = f(X.array, scale_factor)
    NamedArray(S, X.dicts, X.dimnames)
end
