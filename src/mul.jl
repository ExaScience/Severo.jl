import SparseArrays: AbstractSparseMatrixCSC, getcolptr, nonzeros, rowvals, SparseVectorUnion
import LinearAlgebra: mul!, Adjoint, Transpose

#=
mul!(y::AbstractVector, adjA::Adjoint{<:Any,<:AbstractSparseMatrixCSC}, x::SparseVectorUnion, α::Number, β::Number) =
    (A = adjA.parent; _At_or_Ac_mul_B!((a,b) -> adjoint(a) * b, y, A, x, α, β))

function _At_or_Ac_mul_B!(tfun::Function,
                          y::AbstractVector, A::AbstractSparseMatrixCSC, x::SparseVectorUnion,
                          α::Number, β::Number)
    require_one_based_indexing(y, A, x)
    m, n = size(A)
    length(x) == m && length(y) == n || throw(DimensionMismatch())
    n == 0 && return y
    if β != one(β)
        β == zero(β) ? fill!(y, zero(eltype(y))) : rmul!(y, β)
    end
    α == zero(α) && return y

    xnzind = nonzeroinds(x)
    xnzval = nonzeros(x)
    Acolptr = getcolptr(A)
    Arowval = rowvals(A)
    Anzval = nonzeros(A)
    mx = length(xnzind)

    for j = 1:n
        # s <- dot(A[:,j], x)
        s = _spdot(tfun, Acolptr[j], Acolptr[j+1]-1, Arowval, Anzval,
                   1, mx, xnzind, xnzval)
        @inbounds y[j] += s * α
    end
    return y
end
=#

function mul!(y::AbstractVector, A::AbstractSparseMatrixCSC, x::SparseVectorUnion, α::Number, β::Number)
    m, n = size(A)
    length(x) == n && length(y) == m || throw(DimensionMismatch())
    n == 0 && return y
    if β != one(β)
        β == zero(β) ? fill!(y, zero(eltype(y))) : rmul!(y, β)
    end
    α == zero(α) && return y

    nzvA = nonzeros(A)
    rvA = rowvals(A)
    rvx = nonzeroinds(x)
    nzvx = nonzeros(x)
    mx = length(nzvx)

    @inbounds begin
        for jp in 1:mx
            nzx = nzvx[jp]
            j = rvx[jp]
            for kp in nzrange(A, j)
                i = rvA[kp]
                y[i] += α * nzvA[kp] * nzx
            end
        end
    end

    y
end

mul!(C::StridedVecOrMat, A::Transpose{<:Any,<:AbstractSparseMatrixCSC}, B::AbstractSparseMatrixCSC, α::Number, β::Number) = mul!(C, copy(A), B, α, β)
mul!(C::StridedVecOrMat, A::Adjoint{<:Any,<:AbstractSparseMatrixCSC}, B::AbstractSparseMatrixCSC, α::Number, β::Number) = mul!(C, copy(A), B, α, β)

function mul!(C::StridedVecOrMat, A::AbstractSparseMatrixCSC, B::AbstractSparseMatrixCSC, α::Number, β::Number)
    size(A, 1) == size(C, 1) || throw(DimensionMismatch())
    size(A, 2) == size(B, 1) || throw(DimensionMismatch())
    size(B, 2) == size(C, 2) || throw(DimensionMismatch())

    # @inbounds for col in 1:size(B, 2)
    #     mul!(view(C, :, col), A, view(B, :, col), α, β)
    # end

    nzvA = nonzeros(A)
    rvA = rowvals(A)
    nzvB = nonzeros(B)
    rvB = rowvals(B)

    if β != 1
        β != 0 ? rmul!(C, β) : fill!(C, zero(eltype(C)))
    end

    @inbounds begin
        for col in 1:size(B, 2)
            for jp in nzrange(B, col)
                nzB = nzvB[jp]
                j = rvB[jp]
                for kp in nzrange(A, j)
                    row = rvA[kp]
                    C[row, col] += α * nzvA[kp] * nzB
                end
            end
        end
    end

    C
end
