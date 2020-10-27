import Arpack: svds
import LinearAlgebra: svd, Diagonal

struct LinearEmbedding
	coordinates::NamedMatrix
	stdev::NamedVector
	basis::NamedMatrix
end

include("irlba.jl")

function _pca(X; npcs=50, algorithm=:arpack, kw...)
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
	loadings = view(S.V, :, 1:npcs)
	Z, stdev, loadings
end

_pca(X::NamedMatrix; kw...) = _pca(X.array; kw...)
_pca(X::NamedCenteredMatrix; kw...) = _pca(CenteredMatrix(X.A.array, X.mu.array); kw...)

function pca(X::Union{NamedMatrix, NamedCenteredMatrix}; kw...)
	Z, stdev, loadings = _pca(X; kw...)

	k = length(stdev)
	latentnames = map(x -> string("Latent-", x), 1:k)

	rownames, colnames = names(X)
	rowdim, coldim = dimnames(X)

	coordinates = NamedArray(Z, (rownames, latentnames), (rowdim, :latent))
	stdev = NamedArray(stdev, (latentnames,), (:latent,))
	basis = NamedArray(loadings, (colnames, latentnames), (rowdim, :latent))
	LinearEmbedding(coordinates, stdev, basis)
end

function embedding(X; method=:pca, kw...)
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    if method == :pca
        pca(X; kw...)
    else
        error("unknown reduction method: $method")
    end
end

embedding(X::Tuple{<:AbstractMatrix, <:AbstractVector}; kw...) = embedding(CenteredMatrix(first(X), last(X)); kw...)
