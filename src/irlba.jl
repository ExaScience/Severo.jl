# Severo: a software package for analysis and exploration of single-cell RNA-seq datasets.
# Copyright (c) 2021 imec vzw.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, and Additional Terms
# (see below).

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.

import Random: AbstractRNG, default_rng, randn
import LinearAlgebra: mul!, SVD

struct IrlbaData
	rng::AbstractRNG
	A::AbstractMatrix
end

function matmul(yptr::Ptr{Float64}, trans::Cchar, xptr::Ptr{Float64}, temp::Ptr{Float64}, data::IrlbaData)
	A = data.A
	m,n = size(A)

	#T = unsafe_wrap(Array, temp, n)
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

function randv(yptr::Ptr{Float64}, n::BlasInt, data::IrlbaData)
	rng = data.rng
	y = unsafe_wrap(Array, yptr, n)
	y[:] = randn(rng, n)
	nothing
end

function irlba!(rng::AbstractRNG, A::AbstractMatrix, U::Matrix, s::Vector, V::Matrix; init=nothing, tol=1e-5, svtol=tol, maxit=1000, restart=0)
	m,n = size(A)
	nu = length(s)
	m_b = nu + 7

	if m_b < nu
		m_b = nu + 1
	end

	if m_b > min(m,n)
		m_b = min(m,n)
	end

	svtol = min(sqrt(eps()), svtol)

	if init === nothing
		init = randn(rng, n)
	end

	info = ccall(("irlba", libcell), Cint,
		(BlasInt, BlasInt, BlasInt, BlasInt, BlasInt, BlasInt, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid}, Any),
		m, n, nu, m_b, maxit, restart, tol, init, s, U, V,
		@cfunction(randv, Cvoid, (Ptr{Float64}, Cint, Ref{IrlbaData})),
		@cfunction(matmul, Cvoid, (Ptr{Float64}, Cchar, Ptr{Float64}, Ptr{Float64}, Ref{IrlbaData})),
		IrlbaData(rng, A))

	info == 0 || error("convergence failed")

	SVD(U,s,V')
end

function irlba(A::AbstractMatrix, nu::Integer; kw...)
	m,n = size(A)
	V = zeros(n, nu)
	U = zeros(m, nu)
	s = zeros(nu)

	irlba!(default_rng(), A, U, s, V; kw...)
end

function irlba(A::AbstractMatrix, nu::Integer, S::SVD; kw...)
	m,n = size(A)
	V = zeros(n, nu)
	U = zeros(m, nu)
	s = zeros(nu)

	d = length(S.S)
	U[:,1:d] = S.U
	s[1:d] = S.S
	V[:,1:d] = S.V

	irlba!(default_rng(), A, U, s, V; restart=d, kw...)
end
