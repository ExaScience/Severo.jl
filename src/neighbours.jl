# copyright imec - evaluation license - not for distribution

import Random: AbstractRNG, default_rng
import SparseArrays: sparse, nonzeros, droptol!, diag
import Distances: SemiMetric, Euclidean, CosineDist

function ann(rng::AbstractRNG, X::AbstractMatrix{T}, k::Int64, metric::SemiMetric, include_self::Bool=true, ntables::Int64=2*size(X,2)) where {T <: AbstractFloat}
    n,d = size(X)

    nn_index = Matrix{Int32}(undef, n, k)
    distances = Matrix{T}(undef, n, k)

    ann!(rng, nn_index, distances, metric, X, k, include_self, ntables)
end

function ann(rng::AbstractRNG, X::AbstractMatrix{T}, k::Int64; metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where {T <: AbstractFloat}
    ann(rng, X, k, metric, include_self, ntables)
end

function ann!(rng::AbstractRNG, nn_index::Matrix{Int32}, distances::Matrix{Float32}, ::Euclidean, X::StridedMatrix{Float32}, k::Int64, include_self::Bool, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    seed = rand(rng, Int32)
    ccall(("FindNeighboursEuclidean32", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, include_self, seed, nn_index, distances)
    nn_index, distances
end

function ann!(rng::AbstractRNG, nn_index::Matrix{Int32}, distances::Matrix{Float32}, ::CosineDist, X::StridedMatrix{Float32}, k::Int64, include_self::Bool, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    seed = rand(rng, Int32)
    ccall(("FindNeighboursCosine32", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, include_self, seed, nn_index, distances)
    nn_index, distances
end

function ann!(rng::AbstractRNG, nn_index::Matrix{Int32}, distances::Matrix{Float64}, ::Euclidean, X::StridedMatrix{Float64}, k::Int64, include_self::Bool, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    seed = rand(rng, Int32)
    ccall(("FindNeighboursEuclidean64", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, include_self, seed, nn_index, distances)
    nn_index, distances
end

function ann!(rng::AbstractRNG, nn_index::Matrix{Int32}, distances::Matrix{Float64}, ::CosineDist, X::StridedMatrix{Float64}, k::Int64, include_self::Bool, ntables::Int64)
    n,d = size(X)
    @assert stride(X,1) == 1

    seed = rand(rng, Int32)
    ccall(("FindNeighboursCosine64", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, stride(X,2), k, ntables, include_self, seed, nn_index, distances)
    nn_index, distances
end

function _nearest_neighbours(rng::AbstractRNG, X::AbstractMatrix, k::Int64, ::Colon,
        metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2))
    nn_index, distances = ann(rng, X, k, metric, include_self, ntables)
    sparse(vec(nn_index'), repeat(1:size(nn_index,1),inner=k), trues(length(nn_index)))
    #sparse(repeat(1:size(nn_index,1),k), vec(nn_index), trues(length(nn_index)))
end

function _nearest_neighbours(rng::AbstractRNG, X::AbstractMatrix, k::Int64, dims,
        metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2))
    _nearest_neighbours(rng, view(X, :, dims), k, :, metric, include_self, ntables)
end

function _jaccard_index(::Type{T}, nn::SparseMatrixCSC, k, prune::T) where {T <: AbstractFloat}
    snn = convert(SparseMatrixCSC{T}, nn' * nn)
    f(x) = x / (k + (k - x))
    nonzeros(snn) .= f.(nonzeros(snn))
    droptol!(snn, prune)
    snn
end

function _jaccard_index(::Type{T}, nn::SparseMatrixCSC, prune::T) where {T <: AbstractFloat}
    snn = convert(SparseMatrixCSC{T}, nn' * nn)
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
@partial function jaccard_index(nn::NamedArray{T,2} where T; prune::Real=1/15, dtype::Type{R}=Float64) where {R <: AbstractFloat}
    prune = convert(dtype, prune)
    snn = _jaccard_index(dtype, nn.array, prune)
    NamedArray(snn, (nn.dicts[1], nn.dicts[1]), (nn.dimnames[1], nn.dimnames[1]))
end

"""
    jaccard_index(X::NamedArray{T,2}, k::Int64; prune::Real=1/15) where T

Compute a graph with edges defined by the jaccard index. The Jaccard index measures similarity
between nearest neighbour sets, and is defined as the size of the intersection divided by the size
of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `nn`: a nearest neighbour graph
    - `k`: maximum number of neighbours
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
@partial function jaccard_index(nn::NamedArray{T,2} where T, k::Int64; prune::Real=1/15, dtype::Type{R}=Float64) where {R <: AbstractFloat}
    prune = convert(dtype, prune)
    snn = _jaccard_index(dtype, nn.array, k, prune)
    NamedArray(snn, (nn.dicts[1], nn.dicts[1]), (nn.dimnames[1], nn.dimnames[1]))
end

"""
    nearest_neighbours(rng::AbstractRNG, X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell.

**Arguments**:

    - `rng`: random number generator
    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix. k-neighbours are stored as rows for each cell (cols)
"""
function nearest_neighbours(rng::AbstractRNG, X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T
    nn = _nearest_neighbours(rng, X.array, k, dims, metric, include_self, ntables)
    NamedArray(nn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

"""
    nearest_neighbours(X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix. k-neighbours are stored as rows for each cell (cols)
"""
@partial function nearest_neighbours(X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T
    nearest_neighbours(default_rng(), X, k; dims=dims, metric=metric, include_self=include_self, ntables=ntables)
end

"""
    nearest_neighbours(rng::AbstractRNG, em::LinearEmbedding, k::Int64;
         dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on a linear embedding

**Arguments**:

    - `rng`: random number generator
    - `em`: embedding containing the transformed coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix. k-neighbours are stored as rows for each cell (cols)
"""
nearest_neighbours(rng::AbstractRNG, em::LinearEmbedding, k::Int64; dims=:, kw...) = nearest_neighbours(rng, em.coordinates, k; dims=dims, kw...)

"""
    nearest_neighbours(em::LinearEmbedding, k::Int64;
         dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on a linear embedding

**Arguments**:

    - `em`: embedding containing the transformed coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix. k-neighbours are stored as rows for each cell (cols)
"""
nearest_neighbours(em::LinearEmbedding, k::Int64; dims=:, kw...) = nearest_neighbours(default_rng(), em, k; dims=dims, kw...)

"""
    shared_nearest_neighbours(rng::AbstractRNG, X::NamedArray{T,2}, k::Int64; dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `rng`: random number generator
    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
function shared_nearest_neighbours(rng::AbstractRNG, X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2), prune::Real=1/15) where T
    nn = _nearest_neighbours(rng, X.array, k, dims, metric, include_self, ntables)
    prune = convert(T, prune)
    snn = _jaccard_index(T, nn, k, prune)
    NamedArray(snn, (X.dicts[1], X.dicts[1]), (X.dimnames[1], X.dimnames[1]))
end

"""
    shared_nearest_neighbours(rng::AbstractRNG, em::LinearEmbedding, k::Int64; dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2))

Compute a k-nearest neighbours graph based on an embedding of cells and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `rng`: random number generator
    - `em`: embedding containing the transformed coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
shared_nearest_neighbours(rng::AbstractRNG, em::LinearEmbedding, k::Int64; kw...) = shared_nearest_neighbours(rng, em.coordinates, k; kw...)

"""
    shared_nearest_neighbours(X::NamedArray{T,2}, k::Int64;
        dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
@partial shared_nearest_neighbours(X::NamedArray{T,2}, k::Int64; kw...) where T = shared_nearest_neighbours(default_rng(), X, k; kw...)

"""
    shared_nearest_neighbours(em::LinearEmbedding, k::Int64; dims=:, metric::SemiMetric=Euclidean(), include_self::Bool=true, ntables::Int64=2*size(X,2))

Compute a k-nearest neighbours graph based on an embedding of cells and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `em`: embedding containing the transformed coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `dims`: which dimensions to use
    - `include_self`: include the cell in its k-nearest neighbours
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
shared_nearest_neighbours(em::LinearEmbedding, k::Int64; kw...) = shared_nearest_neighbours(default_rng(), em, k; kw...)
