import Random: randn
import LinearAlgebra: mul!, diagm, svd
import CSV
import Test: @test

struct UserData
	A::Matrix{Float64}
end

function matmul(yptr::Ptr{Float64}, trans::Cchar, xptr::Ptr{Float64}, data::UserData)
	m,n = size(data.A)

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

function randv(yptr::Ptr{Float64}, n::Cint, data::UserData)
	y = unsafe_wrap(Array, yptr, n)
	y[:] = randn(n)
	nothing
end

function irlba(A, nu, init=nothing)
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
		@cfunction(randv, Cvoid, (Ptr{Float64}, Cint, Ref{UserData})),
		@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ref{UserData})),
		UserData(A))
	U,s,V
end

function irlba(A, nu, Uold, sold, Vold, init=nothing)
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
		@cfunction(randv, Cvoid, (Ptr{Float64}, Cint, Ref{UserData})),
		@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ref{UserData})),
		UserData(A))
	U,s,V
end

A = Matrix(CSV.read("tests/A1.csv", header=false))

init = [-0.11402883  0.69200073 -0.01087479 -0.09644176  0.06517088 -0.34087181 0.10655937 -0.02089952  0.10182638 -0.21350242 -0.22232905 -0.15070250 -0.16691015  0.28232757  0.34320183 -0.03409794 -0.07800907 -0.07151129 0.01504743  0.02106057]

U,s,V = irlba(A, 2, init)

init = [ 0.075692998 -0.227164504 -0.331333950  0.371658075 -0.186102469 0.057737044 -0.033130097 -0.313616740  0.194434088  0.079816824 0.395964734  0.015153232 -0.096824133  0.165667808 -0.433617654 -0.271884126  0.185128916 -0.151392738 -0.012108516  0.006006806]
U2,s2,V2 = irlba(A, 3, U, s, V, init)
#S = svd(A)
#@test s ≈ S.S[1:nu] atol=1e-5

#S_recon = S.U[:,1:nu] * Diagonal(S.S[1:nu]) * S.V[:,1:nu]'
#A_recon = U*Diagonal(s)*V'
#@test S_recon ≈ A_recon atol=2e-6
