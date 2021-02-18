# copyright imec - evaluation license - not for distribution

import SparseArrays: sparse, nonzeros, droptol!, diag
import Distances: SemiMetric, Euclidean, CosineDist

function ann(X::AbstractMatrix{Float64}, k::Int64, metric::SemiMetric, ntables::Int64=2*size(X,2))
    n,d = size(X)

    nn_index = Matrix{Int32}(undef, n, k)
    distances = Matrix{Float64}(undef, n, k)

    ann!(nn_index, distances, metric, X, k, ntables)
end

function ann!(nn_index::Matrix{Int32}, distances::Matrix{Float64}, ::Euclidean, X::StridedMatrix{Float64}, k::Int64, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    ccall(("FindNeighboursEuclidian", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, nn_index, distances)
    nn_index, distances
end

function ann!(nn_index::Matrix{Int32}, distances::Matrix{Float64}, ::CosineDist, X::StridedMatrix{Float64}, k::Int64, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    ccall(("FindNeighboursCosine", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, nn_index, distances)
    nn_index, distances
end

function nearest_neighbours(X::AbstractMatrix, k::Int64, metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2))
    nn_index, distances = ann(X, k, metric, ntables)
    sparse(repeat(1:size(nn_index,1),k), vec(nn_index), trues(length(nn_index)))
end

function jaccard_index(nn::SparseMatrixCSC, k; prune::Real=1/15)
    snn = convert(SparseMatrixCSC{Float64}, nn * nn')
    f(x) = x / (k + (k - x))
    nonzeros(snn) .= f.(nonzeros(snn))
    droptol!(snn, prune)
    snn
end

function jaccard_index(nn::SparseMatrixCSC; prune::Real=1/15)
    snn = convert(SparseMatrixCSC{Float64}, nn * nn')
    k = diag(snn)

    nz = nonzeros(snn)
    @inbounds for i in 1:size(snn,2)
        x = view(nz, nzrange(snn, i))
        x ./= (k[i] .+ (k[i] .- x))
    end

    droptol!(snn, prune)
    snn
end

"""
    jaccard_index(X::NamedArray{T,2}; prune::Real=1/15) where T

Compute a graph with edges defined by the jaccard index. The Jaccard index measures similarity
between nearest neighbour sets, and is defined as the size of the intersection divided by the size
of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:
    - `nn`: a nearest neighbour graph
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
function jaccard_index(nn::NamedArray{T,2}; prune::Real=1/15) where T
    snn = jaccard_index(nn.array; prune=prune)
    NamedArray(snn, (nn.dicts[1], nn.dicts[1]), (nn.dimnames[1], nn.dimnames[1]))
end

"""
    nearest_neighbours(X::NamedArray{T,2}, k::Int64; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix
"""
function nearest_neighbours(X::NamedArray{T,2}, k::Int64, dims; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2)) where T
    nn = nearest_neighbours(view(X.array,:,dims), k, metric, ntables)
    NamedArray(nn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

function nearest_neighbours(X::NamedArray{T,2}, k::Int64, ::Colon=:;  metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2)) where T
    nn = nearest_neighbours(X.array, k, metric, ntables)
    NamedArray(nn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

"""
    nearest_neighbours(em::LinearEmbedding, k::Int64; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on an embedding

**Arguments**:

    - `em`: embedding containing the transformed coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix
"""
nearest_neighbours(em::LinearEmbedding, k::Int64, dims=:; kw...) = nearest_neighbours(em.coordinates, k, dims; kw...)

"""
    shared_nearest_neighbours(X::NamedArray{T,2}, k::Int64, dims=:; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
function shared_nearest_neighbours(X::NamedArray{T,2}, k::Int64, dims; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2), prune=1/15) where T
    nn = nearest_neighbours(view(X.array, :,dims), k, matric, ntables)
    snn = jaccard_index(nn, k; prune=prune)
    NamedArray(snn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

function shared_nearest_neighbours(X::NamedArray{T,2}, k::Int64, ::Colon=:; metric::SemiMetric=Euclidean(), ntables::Int64=2*size(X,2), prune=1/15) where T
    nn = nearest_neighbours(X.array, k, metric, ntables)
    snn = jaccard_index(nn, k; prune=prune)
    NamedArray(snn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

shared_nearest_neighbours(em::LinearEmbedding, k::Int64, dims=:; kw...) = shared_nearest_neighbours(em.coordinates, k, dims; kw...)
