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

function log_VMR(x::Union{SparseColumnView,SparseVector})
    n = length(x)

    count = n - nnz(x)
    mu = var = zero(eltype(x))

    # nonzeros
    for v in nonzeros(x)
        count += 1
        expm1_v = expm1(v)
        delta = (expm1_v - mu)
        mu += delta / count
        var += delta * (expm1_v - mu)
    end

    log(var/mu)
end

"""
variance to mean ratio (VMR) in non-logspace
"""
function log_VMR(A::SparseMatrixCSC)
    [log_VMR(x) for x in eachcol(A)]
end

function scale_data(A::SparseMatrixCSC; scale_max=10)
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

function show(io::IO, C::CenteredMatrix)
    print(io, "CenteredMatrix(A=$(C.A), mu=$(C.mu))")
end

function show(io::IO, ::MIME"text/plain", C::CenteredMatrix)
    println(io, "CenteredMatrix:")
    ioc = IOContext(io, :compact=>true, :limit=>true)
    println(ioc, "  A = ", C.A)
    print(ioc, "  mu = ", C.mu)
end

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

names(C::NamedCenteredMatrix) = names(C.A)
names(C::NamedCenteredMatrix, d::Integer) = names(C.A, d)
dimnames(C::NamedCenteredMatrix) = dimnames(C.A)
dimnames(C::NamedCenteredMatrix, d::Integer) = dimnames(C.A, d)

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
