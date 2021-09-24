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

function show(io::IO, em::LinearEmbedding; dims::AbstractVector{<:Integer}=1:5, nfeatures::Integer=20)
    loadings = em.basis
    nfeatures = min(nfeatures, size(loadings, 1))
    dims = intersect(dims, axes(loadings, 2))

    for dim in dims
        x = loadings.array[:,dim]
        num = round(Int64, nfeatures / 2)
        positives = partialsortperm(x, 1:num, rev=true)
        negatives = partialsortperm(x, 1:num, rev=false)
        println(io, names(loadings,2)[dim])
        println(io, "Positive: ", names(loadings,1)[positives])
        println(io, "Negative: ", names(loadings,1)[negatives])
    end
end

function show(io::IO, C::CenteredMatrix)
    print(io, "CenteredMatrix(A=$(C.A), mu=$(C.mu))")
end

function show(io::IO, ::MIME"text/plain", C::CenteredMatrix)
    println(io, "CenteredMatrix:")
    ioc = IOContext(io, :compact=>true, :limit=>true)
    println(ioc, "  A = ", C.A)
    print(ioc, "  mu = ", C.mu)
end

