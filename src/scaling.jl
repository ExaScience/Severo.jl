# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, Adjoint, nonzeros, nonzeroinds, nnz, nzrange
import LinearAlgebra: axpy!, mul!

function mean_var(::Type{T}, x::Union{SparseColumnView,SparseVector}) where T
    n = length(x)

    count = n - nnz(x)
    mu = s = zero(T)

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

function mean_std(::Type{T}, x::Union{SparseColumnView,SparseVector}) where T
    mu, var = mean_var(T, x)
    mu, sqrt(var)
end

mean_var(x::Union{SparseColumnView,SparseVector}) = mean_var(Float64, x)
mean_std(x::Union{SparseColumnView,SparseVector}) = mean_std(Float64, x)
mean_var(x::Union{SparseColumnView{T},SparseVector{T}}) where {T<:AbstractFloat} = mean_var(T, x)
mean_std(x::Union{SparseColumnView{T},SparseVector{T}}) where {T<:AbstractFloat} = mean_std(T, x)

function mean_var(::Type{T}, x::Union{SparseColumnView,SparseVector}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) where T
    @assert length(x) == length(lbls)

    mu = zeros(T, nlabels)
    var = zeros(T, nlabels)
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

function mean_var(::Type{T}, x::AbstractVector, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) where T
    @assert length(x) == length(lbls)

    mu = zeros(T, nlabels)
    var = zeros(T, nlabels)
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

function mean_std(::Type{T}, x::Union{SparseColumnView,SparseVector, AbstractVector}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) where T
    mu, var = mean_var(T, x, lbls, nlabels)
    var .= sqrt.(var)
    mu, var
end

mean_var(x::Union{SparseColumnView,SparseVector,AbstractVector}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) = mean_var(Float64, x, lbls, nlabels)
mean_std(x::Union{SparseColumnView,SparseVector,AbstractVector}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) = mean_std(Float64, x, lbls, nlabels)
mean_var(x::Union{SparseColumnView{T},SparseVector{T},AbstractVector{T}}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) where {T<:AbstractFloat} = mean_var(T, x, lbls, nlabels)
mean_std(x::Union{SparseColumnView{T},SparseVector{T},AbstractVector{T}}, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls)) where {T<:AbstractFloat} = mean_std(T, x, lbls, nlabels)

function mean_std(::Type{T}, A::SparseMatrixCSC) where T
    n, d = size(A)
    mu = zeros(T, d)
    std = zeros(T, d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], v = mean_var(T, a)
        std[i] = sqrt(v)
    end

    mu, std
end

function mean_var(::Type{T}, A::SparseMatrixCSC) where T
    n, d = size(A)
    mu = zeros(T, d)
    var = zeros(T, d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], var[i] = mean_var(T, a)
    end

    mu, var
end

mean_std(A::SparseMatrixCSC) = mean_std(Float64, A)
mean_var(A::SparseMatrixCSC) = mean_var(Float64, A)
mean_std(A::SparseMatrixCSC{T}) where {T<:AbstractFloat} = mean_std(T, A)
mean_var(A::SparseMatrixCSC{T}) where {T<:AbstractFloat} = mean_var(T, A)

mean_std(::Type{T}, adjA::Adjoint{<:Any,<:SparseMatrixCSC}) where T = mean_std(T, copy(adjA))
mean_std(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_std(copy(adjA))
mean_var(::Type{T}, adjA::Adjoint{<:Any,<:SparseMatrixCSC}) where T = mean_var(T, copy(adjA))
mean_var(adjA::Adjoint{<:Any,<:SparseMatrixCSC}) = mean_var(copy(adjA))

"""
variance to mean ratio (VMR) in non-logspace
"""
function mean_var(::Type{T}, f::Function, x::Union{SparseColumnView,SparseVector}) where T
    n = length(x)

    count = n - nnz(x)
    mu = var = zero(T)

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
mean_var(f::Function, x::Union{SparseColumnView,SparseVector}) = mean_var(Float64, f, x)

function mean_var(::Type{T}, f::Function, A::SparseMatrixCSC) where T
    n, d = size(A)
    mu = zeros(T, d)
    var = zeros(T, d)

    @inbounds for (i, a) in enumerate(eachcol(A))
        mu[i], var[i] = mean_var(T, f, a)
    end

    mu, var
end
mean_var(f::Function, A::SparseMatrixCSC) = mean_var(Float64, f, A)

function log_VMR(::Type{T}, f::Function, A::SparseMatrixCSC) where T
    mu, var = mean_var(T, f, A)
    log1p.(mu), log.(var ./ mu)
end

log_VMR(::Type{T}, A::SparseMatrixCSC) where T = log_VMR(T, identity, A)
log_VMR(A::SparseMatrixCSC) = log_VMR(Float64, identity, A)
log_VMR(A::SparseMatrixCSC{T}) where {T<:AbstractFloat} = log_VMR(T, identity, A)

function scale_data(A::SparseMatrixCSC, scale_max::R=Inf) where {R <: AbstractFloat}
    n, d = size(A)
    B = similar(A, R)

    # ((x - mu) / std > scale_max)
    mu = zeros(R, d)
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
eltype(::CenteredMatrix{P, T, R}) where {P,T,R} = P
size(S::CenteredMatrix) = size(S.A)
adjoint(S::CenteredMatrix) = Adjoint(S)

_A_mu(S::CenteredMatrix) = S.A, S.mu
_A_mu(S::NamedCenteredMatrix) = S.A.array, S.mu.array

function mul!(C::StridedVector, S::CenteredMatrix, v::StridedVector, α::Number, β::Number)
    A, mu = _A_mu(S)
    mul!(C, A, v, α, β)
    C .-= α * dot(mu, v)
    C
end

function mul!(C::StridedVector, adjS::Adjoint{<:Any, <:CenteredMatrix}, v::StridedVector, α::Number, β::Number)
    S = adjS.parent
    A, mu = _A_mu(S)
    mul!(C, adjoint(A), v,  α, β)
    axpy!(-α * sum(v), mu, C)
end

function mul!(C::StridedMatrix, S::CenteredMatrix, v::StridedMatrix, α::Number, β::Number)
    A, mu = _A_mu(S)
    mul!(C, A, v, α, β)
    C .-= α * mu' * v
    C
end

function mul!(C::StridedMatrix, adjS::Adjoint{<:Any, <:CenteredMatrix}, v::StridedMatrix, α::Number, β::Number)
    S = adjS.parent
    A, mu = _A_mu(S)
    mul!(C, adjoint(A), v,  α, β)
    s = sum(v, dims=1)
    mul!(C, mu, s, α, 1.0)
end

function mul!(C::StridedMatrix, adjS::Adjoint{<:Any, <:CenteredMatrix}, Sr::CenteredMatrix, α::Number, β::Number)
    Sl = adjS.parent
    Al, mul = _A_mu(Sl)
    Ar, mur = _A_mu(Sr)

    # C = A'A - A'M - M'A + M'M

    mul!(C, adjoint(Al), Ar, α, β) # A'A

    q = sum(Ar, dims=1)
    mul!(C, mul, q, -α, 1.0) # - M'A

    q = if Sl !== Sr
        sum(Al, dims=1)
    else
        q
    end

    axpy!(-size(Al,1), mul, q)
    mul!(C, adjoint(q), adjoint(mur), -α, 1.0) # + (M'M - A'M)

    C
end

function convert(::Type{T}, C::CenteredMatrix) where {T<:Array}
    P = eltype(C)
    X = T{P}(C.A)
    X .-= C.mu'
    X
end

function convert(::Type{T}, C::NamedCenteredMatrix) where {T<:Array}
    P = eltype(C)
    X = T{P}(C.A.array)
    X .-= C.mu.array'
    X
end

function convert(::Type{<:NamedArray}, C::NamedCenteredMatrix)
    NamedArray(convert(Matrix, C), C.A.dicts, C.A.dimnames)
end

names(C::NamedCenteredMatrix) = names(C.A)
names(C::NamedCenteredMatrix, d::Integer) = names(C.A, d)
dimnames(C::NamedCenteredMatrix) = dimnames(C.A)
dimnames(C::NamedCenteredMatrix, d::Integer) = dimnames(C.A, d)

"""
    scale_features(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}} ; scale_max=Inf, dtype::Type{<:AbstractFloat})

Scale and center a count/data matrix along the cells such that each feature is standardized

**Arguments**:

    - `X`: the labelled count/data matrix to scale
    - `scale_max`: maximum value of the scaled data

**Return values**:

A centered matrix
"""
function scale_features(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}}; scale_max::Real=Inf, dtype::Type{<:AbstractFloat}=Float64) where T
    scale_max = convert(dtype, scale_max)
    B, mu = scale_data(X.array, scale_max)
    CenteredMatrix(NamedArray(B, X.dicts, X.dimnames), NamedArray(mu, (X.dicts[2],), (X.dimnames[2],)))
end

function scale_features(X::NamedArray{T, 2, SparseMatrixCSC{T, Int64}}; scale_max::Real=Inf, dtype::Type{<:AbstractFloat}=T) where {T <: AbstractFloat}
    scale_max = convert(dtype, scale_max)
    B, mu = scale_data(X.array, scale_max)
    CenteredMatrix(NamedArray(B, X.dicts, X.dimnames), NamedArray(mu, (X.dicts[2],), (X.dimnames[2],)))
end
