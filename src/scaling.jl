# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, Adjoint, nonzeros, nonzeroinds, nnz, nzrange

function mean_var(x::Union{SparseColumnView,SparseVector})
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

    var = s / (n - 1)
    mu, var
end

function mean_std(x::Union{SparseColumnView,SparseVector})
    mu, var = mean_var(x)
    mu, sqrt(var)
end

function mean_var(x::Union{SparseColumnView,SparseVector}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls))
    @assert length(x) == length(lbls)

    mu = zeros(nlabels)
    var = zeros(nlabels)
    n = zeros(Int64, nlabels)
    nz = zeros(Int64, nlabels)

    prev = 1
    @inbounds for i in 1:nnz(x)
        idx = nonzeroinds(x)[i]
        v = nonzeros(x)[i]

        while prev < idx
            k = lbls[prev]
            nz[k] += 1
            prev += 1
        end

        k = lbls[idx]
        n[k] += 1
        delta = (v - mu[k])
        mu[k] += delta / n[k]
        var[k] += delta * (v - mu[k])
        prev = idx + 1
    end

    @inbounds while prev <= length(x)
        k = lbls[prev]
        nz[k] += 1
        prev += 1
    end

    @inbounds for k in 1:nlabels
        delta = mu[k]
        mu[k] *= n[k] / (nz[k] + n[k])
        var[k] += nz[k] * (mu[k]*delta)
        n[k] += nz[k]
    end

    var .= (var ./ (n .- 1))
    mu, var, n
end

function mean_var(x::AbstractVector, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls))
    @assert length(x) == length(lbls)

    mu = zeros(nlabels)
    var = zeros(nlabels)
    n = zeros(Int64, nlabels)

    @inbounds for (v,k) in zip(x, lbls)
        n[k] += 1
        delta = (v - mu[k])
        mu[k] += delta / n[k]
        var[k] += delta * (v - mu[k])
    end

    var ./= (n .- 1)
    mu, var, n
end

function mean_std(x::AbstractVector, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls))
    mu, var = mean_var(x, lbls, nlabels)
    var .= sqrt.(var)
    mu, var
end

function mean_std(A::SparseMatrixCSC)
    n, d = size(A)
    mu = zeros(d)
    std = zeros(d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], v = mean_var(a)
        std[i] = sqrt(v)
    end

    mu, std
end

function mean_var(A::SparseMatrixCSC)
    n, d = size(A)
    mu = zeros(d)
    var = zeros(d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], var[i] = mean_var(a)
    end

    mu, var
end

mean_std(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_std(copy(adjA))
mean_var(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_var(copy(adjA))

"""
variance to mean ratio (VMR) in non-logspace
"""
function mean_var(f::Function, x::Union{SparseColumnView,SparseVector})
    n = length(x)

    count = n - nnz(x)
    mu = var = zero(eltype(x))

    # nonzeros
    for v in nonzeros(x)
        count += 1
        f_v = f(v)
        delta = (f_v - mu)
        mu += delta / count
        var += delta * (f_v - mu)
    end

    var /= (n .- 1)
    mu, var
end

function mean_var(f::Function, A::SparseMatrixCSC)
    n, d = size(A)
    mu = zeros(d)
    var = zeros(d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], var[i] = mean_var(f, a)
    end

    mu, var
end

function log_VMR(f::Function, A::SparseMatrixCSC)
    mu, var = mean_var(f, A)
    log1p.(mu), log.(var ./ mu)
end

log_VMR(A::SparseMatrixCSC) = log_VMR(identity, A)

function scale_data(A::SparseMatrixCSC; scale_max::Float64=10.)
    n, d = size(A)
    B = similar(A, Float64)

    # ((x - mu) / std > scale_max)
    mu = zeros(d)
    @inbounds for i in 1:size(A,2)
        mu[i], std = mean_std(view(A, :, i))
        mu[i] /= std

        scale_max_i = scale_max + mu[i]
        @inbounds for j in nzrange(A, i)
            v = nonzeros(A)[j] / std
            nonzeros(B)[j] = (v > scale_max_i) ? scale_max_i : v
        end
    end

    B, mu
end

struct CenteredMatrix{P, T <: AbstractMatrix, R <: AbstractVector} <: AbstractArray{P,2}
    A::T
    mu::R
end

function CenteredMatrix(A::AbstractMatrix, mu::AbstractVector)
    m,n = size(A)
    @assert n == length(mu)

    T = typeof(A)
    R = typeof(mu)
    P = promote_type(eltype(A),eltype(mu))
    CenteredMatrix{P, T, R}(A, mu)
end

NamedCenteredMatrix{P,T,R} = CenteredMatrix{P,T,R} where {P, T <: NamedArray, R <: AbstractVector}

import Base: eltype, size, adjoint, show, IO, convert
import LinearAlgebra: mul!, dot, axpy!
eltype(S::CenteredMatrix{P, T, R}) where {P,T,R} = P
size(S::CenteredMatrix) = size(S.A)
adjoint(S::CenteredMatrix) = Adjoint(S)

function mul!(C::StridedVector, S::CenteredMatrix, v::StridedVector, α::Number, β::Number)
    mul!(C, S.A, v, α, β)
    C .-= α * dot(S.mu, v)
    C
end

function mul!(C::StridedVector, adjS::Adjoint{<:Any, <:CenteredMatrix}, v::StridedVector, α::Number, β::Number)
    S = adjS.parent
    mul!(C, adjoint(S.A), v,  α, β)
    axpy!(-α * sum(v), S.mu, C)
end

function mul!(C::StridedMatrix, S::CenteredMatrix, v::StridedMatrix, α::Number, β::Number)
    mul!(C, S.A, v, α, β)
    C .-= α * S.mu' * v
    C
end

function mul!(C::StridedMatrix, adjS::Adjoint{<:Any, <:CenteredMatrix}, v::StridedMatrix, α::Number, β::Number)
    S = adjS.parent
    mul!(C, adjoint(S.A), v,  α, β)
    s = sum(v, dims=1)
    mul!(C, S.mu, s, α, 1.0)
end

function convert(::Type{T}, C::CenteredMatrix) where {T<:Array}
    X = (C.A .- C.mu')
    convert(T, X)
end

function convert(::Type{<:NamedArray}, C::NamedCenteredMatrix)
    NamedArray(convert(Matrix, C), C.A.dicts, C.A.dimnames)
end

names(C::NamedCenteredMatrix) = names(C.A)
names(C::NamedCenteredMatrix, d::Integer) = names(C.A, d)
dimnames(C::NamedCenteredMatrix) = dimnames(C.A)
dimnames(C::NamedCenteredMatrix, d::Integer) = dimnames(C.A, d)

"""
    scale_features(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}} ; scale_max=Inf)

Scale and center a count/data matrix along the cells such that each feature is standardized

**Arguments**:

    - `X`: the labelled count/data matrix to scale
    - `scale_max`: maximum value of the scaled data

**Return values**:

A centered matrix
"""
function scale_features(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}}; scale_max::Real=Inf) where T
    scale_max = convert(Float64, scale_max)
    B, mu = scale_data(X.array; scale_max=scale_max)
    CenteredMatrix(NamedArray(B, X.dicts, X.dimnames), NamedArray(mu, (X.dicts[2],), (X.dimnames[2],)))
end
