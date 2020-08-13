import Random: randn
import LinearAlgebra: mul!, diagm, svd
import CSV

A = Matrix(CSV.read("tests/A2.csv", header=false))

m,n = size(A)
nu = 4
m_b = nu + 7

V = zeros(n, nu)
init = [
-0.08399370
 0.21298677
 0.06070336
 0.04599864
-0.08627140
-0.24730293
-0.22866022
 0.19982952
 0.38052247
 0.29374130
 0.05952441
-0.49340295
 0.04274859
 0.04448834
 0.10195045
 0.17584929
-0.10353195
 0.10870144
-0.48196464
-0.05424818
]
U = zeros(m, nu)
s = zeros(nu)

SVtol = min(sqrt(eps()), 1e-6)

function matmul(yptr::Ptr{Float64}, trans::Cchar, xptr::Ptr{Float64}, data::Ptr{Cvoid})
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

function randv(yptr::Ptr{Float64}, n::Cint)
	y = unsafe_wrap(Array, yptr, n)
	y[:] = randn(n)
	nothing
end

restart = 0
tol = 1e-5
maxit = 1000

ccall(("irlba", "./irlba.so"), Cint,
	(Int64, Int64, Int64, Int64, Int64, Int64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
	m, n, nu, m_b, maxit, restart, tol, init, s, U, V,
	@cfunction(randv, Cvoid, (Ptr{Float64}, Cint)),
	@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Cvoid})),
	C_NULL)

S = svd(A)
@test s ≈ S.S[1:nu]

#S_recon = S.U[:,1:nu] * Diagonal(S.S[1:nu]) * S.V[:,1:nu]'
#A_recon = U*Diagonal(s)*V'
#@test S_recon ≈ A_recon atol=2e-6
