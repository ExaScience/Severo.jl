# copyright imec - evaluation license - not for distribution

import Arpack: svds
import LinearAlgebra: svd, Diagonal
import UMAP

struct LinearEmbedding
    parent::Union{NamedMatrix, NamedCenteredMatrix}
    coordinates::NamedMatrix
    stdev::NamedVector
    basis::NamedMatrix
end

include("irlba.jl")

function _pca(X, npcs::Int64; algorithm=:irlba, kw...)
    m,n = size(X)
    npcs = min(min(m,n), npcs)

    if (npcs > 0.5 * min(m, n)) && (algorithm == :arpack || algorithm == :irlba)
        @warn "Computing too large a percentage of principal components, using standard svd instead"
        algorithm = :svd
    end

    S = if algorithm == :arpack
        S, nconv, niter, nmult, resid = svds(X; nsv=npcs, kw...)
        S
    elseif algorithm == :irlba
        irlba(X, npcs; kw...)
    else
        Q = convert(Matrix, X)
        svd(Q; kw...)
    end

    Z = view(S.U, :, 1:npcs) * Diagonal(view(S.S, 1:npcs))
    stdev = view(S.S, 1:npcs) ./ sqrt(max(1, size(X,1) - 1))
    loadings = if npcs != size(S.V, 2)
        view(S.V, :, 1:npcs)
    else
        S.V
    end

    Z, stdev, loadings
end

_pca(X::NamedMatrix, npcs::Int64; kw...) = _pca(X.array, npcs; kw...)
_pca(X::NamedCenteredMatrix, npcs::Int64; kw...) = _pca(CenteredMatrix(X.A.array, X.mu.array), npcs; kw...)

function pca(X::Union{NamedMatrix, NamedCenteredMatrix}, npcs::Int64; kw...)
    Z, stdev, loadings = _pca(X, npcs; kw...)

    k = length(stdev)
    latentnames = map(x -> string("PC-", x), 1:k)

    rownames, colnames = names(X)
    rowdim, coldim = dimnames(X)

    coordinates = NamedArray(Z, (rownames, latentnames), (rowdim, :latent))
    stdev = NamedArray(stdev, (latentnames,), (:latent,))
    basis = NamedArray(loadings, (colnames, latentnames), (rowdim, :latent))
    LinearEmbedding(X, coordinates, stdev, basis)
end

function UMAP.knn_search(X::AbstractMatrix, k, ::Val{:ann})
    @time knns, dists = ann(X', k, CosineDist(), false)
    knns', dists'
end

"""
    umap(em::LinearEmbedding, ncomponents::Int64=2; dims=:, metric=:cosine, nneighbours::Int=30, min_dist::Real=.3, nepochs::Int=300, kw...) where T

Performs a Uniform Manifold Approximation and Projection (UMAP) dimensional reduction on the coordinates in the linear embedding.

For a more in depth discussion of the mathematics underlying UMAP, see the ArXiv paper: [https://arxiv.org/abs/1802.03426]

**Arguments**:

    - `em`: embedding containing the transformed coordinates for each cell
    - `ncomponents`: the dimensionality of the embedding
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `nneighbours`: the number of neighboring points used in local approximations of manifold structure.
    - `min_dist`: controls how tightly the embedding is allowed compress points together.
    - `nepochs`: number of training epochs to be used while optimizing the low dimensional embedding
    - `kw`: additional parameters for the umap algorithm. See [`UMAP.umap`](@ref)

**Return values**:

A low-dimensional embedding of the cells
"""
umap(em::LinearEmbedding, ncomponents::Int64=2,; dims=:, kw...) = _umap(em, ncomponents, dims; kw...)
_umap(em::LinearEmbedding, ncomponents::Int64, ::Colon; kw...) = umap(em.coordinates, ncomponents; kw...)
_umap(em::LinearEmbedding, ncomponents::Int64, dims; kw...) = umap(view(em.coordinates,:,dims), ncomponents; kw...)

function umap(X::AbstractMatrix, ncomponents::Int64=2; metric=:cosine, nneighbours::Integer=30, min_dist::Real=0.3, nepochs::Integer=300, kw...)
    metric = if metric == :cosine
        UMAP.CosineDist()
    elseif metric == :euclidian
        UMAP.Euclidian()
    else
        metric
    end

    UMAP.umap(X', ncomponents; metric=metric, n_neighbors=nneighbours, min_dist=min_dist, n_epochs=nepochs, kw...)'
end

"""
    umap(X::NamedMatrix, ncomponents::Int64=2; dims=:, metric=:cosine, nneighbours::Int=30, min_dist::Real=.3, nepochs::Int=300, kw...) where T

Performs a Uniform Manifold Approximation and Projection (UMAP) dimensional reduction on the coordinates.

For a more in depth discussion of the mathematics underlying UMAP, see the ArXiv paper: [https://arxiv.org/abs/1802.03426]

**Arguments**:

    - `X`: a labelled matrix with coordinates for each cell
    - `ncomponents`: the dimensionality of the embedding
    - `dims`: which dimensions to use
    - `metric`: distance metric to use
    - `nneighbours`: the number of neighboring points used in local approximations of manifold structure.
    - `min_dist`: controls how tightly the embedding is allowed compress points together.
    - `nepochs`: number of training epochs to be used while optimizing the low dimensional embedding
    - `kw`: additional parameters for the umap algorithm. See [`UMAP.umap`](@ref)

**Return values**:

A low-dimensional embedding of the cells
"""
function umap(X::NamedMatrix, ncomponents::Int64=2; kw...)
    coords = umap(X.array, ncomponents; kw...)

    rownames = names(X, 1)
    rowdim = dimnames(X, 1)
    latentnames = map(x -> string("UMAP-", x), 1:ncomponents)
    NamedArray(coords, (rownames, latentnames), (rowdim, :latent))
end

function embedding(X, ncomponents::Int64; method=:pca, kw...)
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    if method == :pca
        pca(X, ncomponents; kw...)
    elseif method == :umap
        umap(X, ncomponents; kw...)
    else
        error("unknown reduction method: $method")
    end
end
