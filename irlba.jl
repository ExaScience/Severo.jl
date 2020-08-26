import Random: randn
import LinearAlgebra: mul!, diagm, svd, dot
import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange, sprand

abstract type MatMulData end

struct MatMulDataNoCenter{T} <: MatMulData
	A::T
end

function matmul(yptr::Ptr{Float64}, trans::Cchar, xptr::Ptr{Float64}, temp::Ptr{Float64}, data::MatMulDataNoCenter)
	A = data.A
	m,n = size(A)

	if trans == 84
		x = unsafe_wrap(Array, xptr, m)
		y = unsafe_wrap(Array, yptr, n)
		mul!(y, A', x)
	else
		x = unsafe_wrap(Array, xptr, n)
		y = unsafe_wrap(Array, yptr, m)
		mul!(y, A, x)
	end
	nothing
end

struct MatMulDataCenterScale{T} <: MatMulData
	A::T
	center::Vector{Float64}
	scale::Vector{Float64}
end

function matmul(yptr::Ptr{Float64}, trans::Cchar, xptr::Ptr{Float64}, temp::Ptr{Float64}, data::MatMulDataCenterScale)
	A = data.A
	m,n = size(A)

	if trans == 84
		x = unsafe_wrap(Array, xptr, m)
		y = unsafe_wrap(Array, yptr, n)
		mul!(y, A', x)
		y ./= data.scale

		beta = sum(x)
		y .-= beta .* data.center ./ data.scale
	else
		x = unsafe_wrap(Array, xptr, n)
		y = unsafe_wrap(Array, yptr, m)
		T = unsafe_wrap(Array, temp, n)
		T .= x ./ data.scale
		mul!(y, A, T)
		y .-= dot(data.center, T)
	end
	nothing
end

function MatMulData(A::T, center=nothing, scale=nothing) where T
	if center !== nothing
		MatMulDataCenterScale(A, center, scale)
	else
		MatMulDataNoCenter(A)
	end
end

function randv(yptr::Ptr{Float64}, n::Cint, data::MatMulData)
	y = unsafe_wrap(Array, yptr, n)
	y[:] = randn(n)
	nothing
end

function irlba(A, nu; center=nothing, scale=nothing, init=nothing)
	m,n = size(A)
	m_b = nu + 7

	if m_b < nu
		m_b = nu + 1
	end

	if m_b > min(m,n)
		m_b = min(m,n)
	end

	V = zeros(n, nu)
	if init === nothing
		init = randn(n)
	end

	U = zeros(m, nu)
	s = zeros(nu)

	SVtol = min(sqrt(eps()), 1e-6)
	restart = 0
	tol = 1e-5
	maxit = 1000

	ccall(("irlba", "./irlba.so"), Cint,
		(Int64, Int64, Int64, Int64, Int64, Int64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Any),
		m, n, nu, m_b, maxit, restart, tol, init, s, U, V,
		@cfunction(randv, Cvoid, (Ptr{Float64}, Cint, Ref{MatMulData})),
		@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Float64}, Ref{MatMulData})),
		MatMulData(A, center, scale))
	U,s,V
end

function irlba(A, nu, Uold, sold, Vold; center=nothing, scale=nothing, init=nothing)
	m,n = size(A)
	m_b = nu + 7

	if m_b < nu
		m_b = nu + 1
	end

	if m_b > min(m,n)
		m_b = min(m,n)
	end

	V = zeros(n, nu)
	if init === nothing
		init = randn(n)
	end

	U = zeros(m, nu)
	s = zeros(nu)

	d = length(sold)
	U[:,1:d] = Uold
	s[1:d] = sold
	V[:,1:d] = Vold

	SVtol = min(sqrt(eps()), 1e-6)
	restart = d
	tol = 1e-5
	maxit = 1000

	ccall(("irlba", "./irlba.so"), Cint,
		(Int64, Int64, Int64, Int64, Int64, Int64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Any),
		m, n, nu, m_b, maxit, restart, tol, init, s, U, V,
		@cfunction(randv, Cvoid, (Ptr{Float64}, Cint, Ref{MatMulData})),
		@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Float64}, Ref{MatMulData})),
		MatMulData(A, center, scale))
	U,s,V
end
