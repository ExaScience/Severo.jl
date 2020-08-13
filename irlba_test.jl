import Random: randn
import LinearAlgebra: mul!, diagm
import CSV

A = Matrix(CSV.read("A.csv", header=false))

m,n = size(A)
nu = 4
m_b = nu + 7

V = zeros(n, nu)
V[:,1] = [
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
	(Int64, Int64, Int64, Int64, Int64, Int64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
	m, n, nu, m_b, maxit, restart, tol, s, U, V,
	@cfunction(randv, Cvoid, (Ptr{Float64}, Cint)),
	@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Cvoid})),
	C_NULL)

Vout = [
	-0.70265   0.10696  -0.41588   0.23364
	-0.61396  -0.31576  -0.02635  -0.23876
	-0.13002  -0.33378   0.36526  -0.69772
	-0.22419  -0.22507   0.74020   0.58828
	 0.24933  -0.85252  -0.38084   0.23561
]

Wout = [
	-0.08372   0.05147   0.75445   0.58708
	-0.11191   0.39449   0.40767  -0.59833
	 0.61119  -0.20512   0.19593   0.05857
	-0.02407  -0.49257   0.38140  -0.39264
	-0.44195   0.29208  -0.08637   0.24070
	 0.13643  -0.17449  -0.25603   0.09733
	 0.43744   0.08534  -0.04920   0.18957
	 0.03355   0.11618  -0.05446  -0.05063
	 0.44369   0.58077   0.04388  -0.13140
	-0.05530  -0.28843   0.02081  -0.12871
]

Bout = [
	 1.84751   1.32549   0.00000   0.00000
	 0.00000   2.57344   2.01139   0.00000
	 0.00000   0.00000   2.37540   1.21943
	 0.00000   0.00000   0.00000   1.94601
]

Fout = [
	 0.346946
	-0.457926
	 0.336514
	 0.047937
	 0.068736
]
