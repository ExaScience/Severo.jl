import Random: randn
import LinearAlgebra: mul!, diagm, svd, dot
import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange, sprand
import CSV
import Test: @test

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
		y .-= beta .* center ./ scale
	else
		x = unsafe_wrap(Array, xptr, n)
		y = unsafe_wrap(Array, yptr, m)
		temp .= data.scale .* x
		mul!(y, A, temp)
		y .-= dot(data.center, temp)
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

function mean_var(x::Union{SparseColumnView, SparseVector})
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

	std = sqrt(s / (n-1))
	mu, std
 end

 function mean_var(A::SparseMatrixCSC{T}) where T
	n = size(A, 2)
	mu = zeros(T, n)
	std = zeros(T, n)

	for (i,x) in enumerate(eachcol(A))
		mu[i], std[i] = mean_var(x)
	end

	mu, std
 end

#A = Matrix(CSV.read("tests/A1.csv", header=false))

#U,s,V = irlba(A, 2)

#U2,s2,V2 = irlba(A, 3, U, s, V, init)
#S = svd(A)
#@test s ≈ S.S[1:nu] atol=1e-5

#S_recon = S.U[:,1:nu] * Diagonal(S.S[1:nu]) * S.V[:,1:nu]'
#A_recon = U*Diagonal(s)*V'
#@test S_recon ≈ A_recon atol=2e-6

A = sprand(10000, 500, 0.02)
@time U,s,V = irlba(A, 10)
@time S = svd(Matrix(A))
@time UU,ss,VV = irlba(A'A, 10)
