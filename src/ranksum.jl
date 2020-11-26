import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange
import Distributions: cdf, ccdf, Normal

function ranksum(x::Union{SparseColumnView, SparseVector}, labels::BitVector)
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

function ranksum(x::AbstractArray, labels::BitVector)
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

function ranksum(x::AbstractArray, labels::Vector{<:Integer}, nlabels=length(unique(labels)))
    n = length(x)
    J = sortperm(x)

    ri = zeros(Float64, nlabels)
    ni = count_map(labels, nlabels)
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

function ranksumtest(x::AbstractArray, labels::BitVector)
    n = length(x)

    r1, n1, ties_adjust = ranksum(x, labels)

    n2 = n-n1
    U = n1*n2 + n1*(n1+1)/2 - r1

    ties_adjust /= n*(n-1)
    sigma = sqrt((n1 * n2 / 12) * ((n + 1) - ties_adjust))

    mu = n1*n2/2.0

    z = (U - 0.5*sign(U - mu) - mu) / sigma
    twosided = 2 * ccdf(Normal(), abs(z))

    z, twosided
end

function ranksumtest(x::AbstractArray, labels::Vector{<:Integer}, nlabels=length(unique(labels)))
    n = length(x)

    ri, ni, ties_adjust = ranksum(x, labels, nlabels)

    nothers = n .- ni
    U = ni .* nothers + ni .* (ni .+ 1) ./ 2 .- ri

    ties_adjust /= n*(n-1)
    sigma = sqrt.((ni .* nothers ./ 12) .* ((n + 1) - ties_adjust))

    mu = ni.*nothers./2

    zlowertail = (U.+0.5.-mu)./sigma
    zuppertail = (U.-0.5.-mu)./sigma

    z = (U .- 0.5*sign.(U .- mu) .- mu) ./ sigma
    twosided = 2 * ccdf.(Normal(), abs.(z))

    z, twosided
end
