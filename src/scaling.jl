import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, Adjoint, nonzeros, nonzeroinds, nnz, nzrange

function mean_std(x::Union{SparseColumnView,SparseVector})
    n = length(x)

    count = n - nnz(x)
    mu = s = zero(eltype(x))

   # nonzeros
    for v in nonzeros(x)
        count += 1
        delta = (v - mu)
        mu += delta / count
        s += delta * (v - mu)
    end

    std = sqrt(s / (n - 1))
    mu, std
end

function mean_std(A::SparseMatrixCSC)
    n, d = size(A)
    mu = zeros(d)
    std = zeros(d)

    for (i, a) in enumerate(eachcol(A))
        mu[i], std[i] = mean_std(a)
    end

    mu, std
end

function mean_var(A::SparseMatrixCSC)
    n, d = size(A)
    mu = zeros(d)
    var = zeros(d)

    for (i, a) in enumerate(eachcol(A))
        mu[i], std = mean_std(a)
        var[i] = std^2
    end

    mu, var
end

mean_std(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_std(copy(adjA))
mean_var(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_var(copy(adjA))

function scale_data(A::SparseMatrixCSC; scale_max=Inf)
    n, d = size(A)
    B = similar(A, Float64)

    # ((x - mu) / std > scale_max)
    mu = zeros(d)
    @inbounds for i in 1:size(A,2)
        mu[i], std = mean_std(A[:,i])
        mu[i] /= std

        scale_max_i = scale_max + mu[i]
        @inbounds for j in nzrange(A, i)
            v = nonzeros(A)[j] / std
            nonzeros(B)[j] = (v > scale_max_i) ? scale_max_i : v
        end
    end

    B, mu
end

"""
    scale(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}} ; scale_max=Inf)

Scale and center a count/data matrix along the cells such that each feature is standardized

**Arguments**:

    - `X`: the labelled count/data matrix to scale
    - `scale_max`: maximum value of the scaled data

**Return values**:

A labelled scaled matrix and its mean
"""
function scale(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}} ; scale_max=Inf) where T
    B, mu = scale_data(X.array; scale_max=scale_max)
    NamedArray(B, X.dicts, X.dimnames), NamedArray(mu, (X.dicts[2],), (X.dimnames[2],))
end
