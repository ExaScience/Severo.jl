import SparseArrays: SparseMatrixCSC, SparseVector, nonzeros, nonzeroinds, nnz
import Distributions: cdf, ccdf, Normal

function rank(x::AbstractArray)
	n = length(x)
	J = sortperm(x)

	rk = zeros(size(J))
	i = 1
	while i <= n
		j = i
		while (j < n) && x[J[j]] == x[J[j+1]]
			j += 1
		end

		for k in i:j
			rk[J[k]] = (i + j) / 2.
		end

		i = j + 1
	end

	rk
end

function ranksum(x::SparseVector, labels::Vector{Bool})
	n = nnz(x)
	nz = nonzeros(x)

	J = sortperm(nz)
	ii = nonzeroinds(x)[J]
	nzeros = length(x) - n

	r1 = 0.0
	n1 = sum(labels)
	ties_adjust = 0.0

	# sum(rk .* labels) = sum(rk_lt0 .* labels_lt0) + sum(rk_0 .* labels_0) + sum(rk_gt0 .* labels_gt0)
	#   = sum(rk_lt0 .* labels_lt0) + rk_0 * sum(labels_0) + sum((rk_gt0 - rk_0) .* labels_gt0) + rk_0 * sum(labels_gt0)
	#   = sum(rk_lt0 .* labels_lt0) + rk_0 * (sum(labels) - sum(labels_lt0)  + sum((rk_gt0 - rk_0) .* labels_gt0)

	nprezero = 0
	i = 1
	while i <= n && nz[J[i]] < zero(eltype(x))
		j = i
		while (j < n) && nz[J[j]] == nz[J[j+1]]
			j += 1
		end

		rk = (i + j) / 2 # (average) rank for i:j
		r1 += sum(labels[ii[i:j]] .* rk)
		nprezero += sum(labels[ii[i:j]])

		ti = length(i:j)
		ties_adjust += (ti^3 - ti)

		i = j + 1
	end

	zerosrank = i
	rk_0 = zerosrank + (nzeros - 1) / 2

	while i <= n
		j = i
		while (j < n) && nz[J[j]] == nz[J[j+1]]
			j += 1
		end

		rk = (i + j) / 2 - (zerosrank - (nzeros+1)/2)
		r1 += sum(labels[ii[i:j]] .* rk)

		ti = length(i:j)
		ties_adjust += (ti^3 - ti)

		i = j + 1
	end

	r1 += rk_0 * (n1 - nprezero)
	ties_adjust += (nzeros^3 - nzeros)

	r1, n1, ties_adjust
end

function ranksum(x::AbstractArray, labels::Vector{Bool})
	n = length(x)
	J = sortperm(x)

	r1 = 0.0
	n1 = sum(labels)
	ties_adjust = 0.0

	i = 1
	while i <= n
		j = i
		while (j < n) && x[J[j]] == x[J[j+1]]
			j += 1
		end

		rk = (i + j) / 2 # (average) rank for i:j
		r1 += sum(labels[J[i:j]] .* rk)

		ti = length(i:j)
		ties_adjust += (ti^3 - ti)

		i = j + 1
	end

	r1, n1, ties_adjust
end

function counts(x::Vector{T}, nlevels=length(unique(x))) where {T <: Signed}
	n = zeros(Int64, nlevels)
	for i in x
		n[i] += 1
	end
	n
end

function ranksum(x::AbstractArray, labels::Vector{T}, nlabels=length(unique(labels))) where {T <: Signed}
	n = length(x)
	J = sortperm(x)

	ri = zeros(Float64, nlabels)
	ni = counts(labels, nlabels)
	ties_adjust = 0.0

	i = 1
	while i <= n
		j = i
		while (j < n) && x[J[j]] == x[J[j+1]]
			j += 1
		end

		rk = (i + j) / 2 # (average) rank for i:j
		for lbl in labels[J[i:j]]
			ri[lbl] += rk
		end

		ti = length(i:j)
		ties_adjust += (ti^3 - ti)

		i = j + 1
	end

	ri, ni, ties_adjust
end

function ranksumtest(A::SparseMatrixCSC, labels::Vector{Bool})
	n = size(A,1)
	@assert n == length(labels)

	[ranksumtest(A[:,c], labels) for c in 1:size(A,2)]
end

function ranksumtest(x::AbstractArray, labels::Vector{Bool})
	n = length(x)

	r1, n1, ties_adjust = ranksum(x, labels)

	n2 = n-n1
	U = n1*n2 + n1*(n1+1)/2 - r1

	ties_adjust /= n*(n-1)
	sigma = sqrt((n1 * n2 / 12) * ((n + 1) - ties_adjust))

	mu = n1*n2/2.0

	zlowertail = (U+0.5-mu)/sigma
	z = (U - 0.5*sign(U) - mu) / sigma
	zuppertail = (U-0.5-mu)/sigma

	less, greater = cdf(Normal(), zlowertail), ccdf(Normal(), zuppertail)
	twosided = 2 * min(less, greater)

	twosided
end

function ranksumtest(x::AbstractArray, labels::Vector{T}, nlabels=length(unique(labels))) where {T <: Signed}
	n = length(x)

	ri, ni, ties_adjust = ranksum(x, labels, nlabels)

	nothers = n .- ni
	U = ni .* nothers + ni .* (ni .+ 1) ./ 2 .- ri

	ties_adjust /= n*(n-1)
	sigma = sqrt.((ni .* nothers ./ 12) .* ((n + 1) - ties_adjust))

	mu = ni.*nothers./2

	zlowertail = (U.+0.5.-mu)./sigma
	zuppertail = (U.-0.5.-mu)./sigma

	less, greater = cdf.(Normal(), zlowertail), ccdf.(Normal(), zuppertail)
	twosided = 2 .* min.(less, greater)

	twosided
end

function ranksumtest(A::SparseMatrixCSC, labels::Vector{T}, nlabels=length(unique(labels))) where {T <: Signed}
	n = size(A,1)
	@assert n == length(labels)

	[ranksumtest(Vector(A[:,c]), labels, nlabels) for c in 1:size(A,2)]
end

statistics = [1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30, 0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29]
labels = Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
pval = 0.13291945818531886
