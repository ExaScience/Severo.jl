import SparseArrays: sparse, nonzeros, droptol!, diag

function ann(X, k; ntables=2*size(X,2))
    n,d = size(X)

    nn_index = Matrix{Int32}(undef, n, k)
    distances = Matrix{Float64}(undef, n, k)

    ccall(("FindNeighbours", libcell), Cvoid,
        (Ptr{Float64}, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
        X, n, d, k, ntables, nn_index, distances)
    nn_index, distances
end

function nearest_neigbours(X::AbstractMatrix, k::Int64; ntables::Int64=2*size(X,2))
    nn_index, distances = ann(X, k, ntables=ntables)
    sparse(repeat(1:size(nn_index,1),k), vec(nn_index), trues(length(nn_index)))
end

function jaccard_index(nn::SparseMatrixCSC, k; prune=1/15)
    snn = convert(SparseMatrixCSC{Float64}, nn * nn')
    f(x) = x / (k + (k - x))
    nonzeros(snn) .= f.(nonzeros(snn))
    droptol!(snn, prune)
    snn
end

function jaccard_index(nn::SparseMatrixCSC; prune=1/15)
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
    nearest_neigbours(X::NamedArray{T,2}, k::Int64; ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)

**Return values**:

A k-nearest neighbours graph represented by a sparse matrix
"""
function nearest_neigbours(X::NamedArray{T,2}, k::Int64; ntables::Int64=2*size(X,2)) where T
    nn = nearest_neigbours(X.array, k; ntables=ntables)
    NamedArray(nn, (X.dicts[1], X.dicts[1]), (:cells, :cells))
end

nearest_neigbours(em::LinearEmbedding, k::Int64; kw...) = nearest_neigbours(em.coordinates, k; kw...)

"""
    shared_nearest_neigbours(X::NamedArray{T,2}, k::Int64; ntables::Int64=2*size(X,2)) where T

Compute a k-nearest neighbours graph based on coordinates for each cell and its Jaccard index.\\
The Jaccard index measures similarity between nearest neighbour sets, and is defined as
the size of the intersection divided by the size of the union. "0" indicating no overlap and "1" indicating full overlap.

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `k`: number of nearest neighbours to find
    - `ntables`: number of tables to use in knn algorithm: controls the precision (higher is more accurate)
    - `prune`: cutoff for the Jaccard index, edges with values below this cutoff are removed from the resulting graph

**Return values**:

A shared nearest neighbours graph represented by a sparse matrix. Weights of the edges indicate similarity of
the neighbourhoods of the cells as computed with the Jaccard index.
"""
function shared_nearest_neigbours(X::NamedArray{T,2}, k::Int64; ntables::Int64=2*size(X,2), prune=1/15) where T
    nn = nearest_neigbours(X.array, k; ntables=ntables)
    snn = jaccard_index(nn, k; prune=prune)
    NamedArray(snn, (X.dicts[1], X.dicts[1]), (:cells, :cells))
end

shared_nearest_neigbours(em::LinearEmbedding, k::Int64; kw...) = shared_nearest_neigbours(em.coordinates, k; kw...)
