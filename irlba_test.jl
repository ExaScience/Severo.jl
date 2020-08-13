import Random: randn
import LinearAlgebra: mul!, diagm
import CSV

A = Matrix(CSV.read("A.csv", header=false))

m,n = size(A)
nu = 4
m_b = nu + 7

VT = zeros(n, nu)
VT[:,1] = [
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
	m, n, nu, m_b, maxit, restart, tol, s, U, VT,
	@cfunction(randv, Cvoid, (Ptr{Float64}, Cint)),
	@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Cvoid})),
	C_NULL)

VTout = [
	-0.2545375 -0.07467720 -0.346707051 -0.32328830
	-0.2077932  0.16399732 -0.091285509  0.42002493
	-0.2589141 -0.11836489  0.059818607  0.30175572
	-0.1534967  0.09946021 -0.147369671 -0.02521854
	-0.2161617  0.28072620 -0.164620343 -0.14767862
	-0.2530343 -0.29908925  0.038007320  0.30773296
	-0.2371615 -0.02194428 -0.012802755  0.05997200
	-0.1693965  0.17450458 -0.134799368  0.29348798
	-0.2208695  0.14312412  0.012781744 -0.17417566
	-0.2105408 -0.32743596  0.128399353 -0.18794706
	-0.1716069 -0.06786656  0.135366567  0.10332735
	-0.2132134  0.16620796  0.005567743 -0.28333757
	-0.2059571  0.40695363  0.379648759 -0.14270936
	-0.2371838 -0.30787960  0.315891322 -0.28371695
	-0.2176964  0.10316317 -0.375943930  0.16543498
	-0.2411823 -0.32639392  0.071684741  0.21354967
	-0.2395303  0.26690038  0.516584192  0.05414506
	-0.2383706 -0.30132329 -0.199179982 -0.25973748
	-0.2605490  0.21153398 -0.253082862 -0.11391703
	-0.2267290  0.01208149  0.057490410  0.08548152
]

Uout = [
	-0.1930165  0.18231318 -0.28774161  0.202746917
	-0.2263395 -0.11706939 -0.20346207 -0.155162818
	-0.2552328  0.17003377  0.17067742  0.287137432
	-0.2331016 -0.37977598 -0.03424417  0.033210403
	-0.2791854  0.10892746  0.02116594 -0.253240453
	-0.1885287 -0.09726451 -0.02578441  0.698484735
	-0.2680791 -0.14919356 -0.43344269 -0.035343538
	-0.2240754  0.10327027  0.09161648 -0.301769966
	-0.2030925 -0.11751559 -0.02921655 -0.027701834
	-0.1949527 -0.23424844  0.16641441 -0.135215554
	-0.2013967  0.40850072 -0.03887190 -0.314315868
	-0.2025152  0.46620800  0.16075185  0.074105526
	-0.2271294 -0.13337367  0.44515600  0.095385539
	-0.2019736 -0.03818759 -0.24856436 -0.049840011
	-0.2107654 -0.11246333  0.43831198 -0.025275227
	-0.2040979 -0.18749306 -0.22167621  0.018682548
	-0.2274341 -0.28922304 -0.01790059 -0.149711163
	-0.2020057  0.14834484 -0.06512042 -0.064199552
	-0.2605427 -0.02054547  0.25050037 -0.007472995
	-0.2382664  0.31564950 -0.15397597  0.210816335
]

Sout = [10.313760  2.218513  1.987749  1.932419]