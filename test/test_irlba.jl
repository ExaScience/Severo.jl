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

using Severo, Test
import SparseArrays: sprand
import LinearAlgebra: svd, SVD, Diagonal, norm
import Statistics: mean
import Distributions: rand, Poisson

relative_error(X::AbstractMatrix, S::SVD, k::Integer=length(S.S)) = norm(X - S.U[:,1:k] * Diagonal(S.S[1:k]) * S.Vt[1:k,:])
svd(C::Severo.CenteredMatrix) = svd(convert(Matrix, C))

@testset "irlba" begin
	@testset "dense square" begin
		X = randn(50, 50)
		S = Severo.irlba(X, 20; tol=1e-5)
		SS = svd(X)
		@test S.S ≈ SS.S[1:20]
		@test norm(X'*S.U - S.V*Diagonal(S.S))/norm(X) < 1e-5
		@test relative_error(X, S) ≈ relative_error(X, SS, 20)
	end

	@testset "dense non-square" begin
		X = randn(100, 50)
		S = Severo.irlba(X, 20; tol=1e-5)
		SS = svd(X)
		@test S.S ≈ SS.S[1:20]
		@test norm(X'*S.U - S.V*Diagonal(S.S))/norm(X) < 1e-5
		@test relative_error(X, S) ≈ relative_error(X, SS, 20)
	end

	@testset "sparse" begin
		X = sprand(2000, 400, .1)
		S = Severo.irlba(X, 2; tol=1e-9)
		SS = svd(Matrix(X))
		@test S.S ≈ SS.S[1:2]
		@test relative_error(X, S) ≈ relative_error(X, SS, 2)
		@test norm(X'*S.U - S.V*Diagonal(S.S))/norm(X) < 1e-9
	end

	@testset "restart" begin
		X = randn(20, 10)
		SS = svd(X)

		S = Severo.irlba(X, 2; tol=1e-5)
		@test S.S ≈ SS.S[1:2]
		@test relative_error(X, S) ≈ relative_error(X, SS, 2)

		S = Severo.irlba(X, 3, S; tol=1e-5)
		@test_broken S.S[3] ≈ SS.S[3]
		@test_broken relative_error(X, S) ≈ relative_error(X, SS, 3)
	end

	@testset "centered dense matrix" begin
		X = randn(20, 10)
		mu = vec(mean(X, dims=1))

		C = Severo.CenteredMatrix(X, mu)
		Q = convert(Matrix, C)
		SS = svd(Q)

		S = Severo.irlba(C, 3)
		@test S.S ≈ SS.S[1:3]
		@test relative_error(Q, S) ≈ relative_error(Q, SS, 3)
		@test norm(Q'*S.U - S.V*Diagonal(S.S))/norm(Q) < 1e-9
	end

	@testset "centered sparse matrix" begin
		X = sprand(2000, 400, .1)
		mu = vec(mean(X, dims=1))

		C = Severo.CenteredMatrix(X, mu)
		Q = convert(Matrix, C)
		SS = svd(Q)

		S = Severo.irlba(C, 2, tol=1e-9)
		@test S.S ≈ SS.S[1:2]
		@test relative_error(Q, S) ≈ relative_error(Q, SS, 2)
		@test norm(Q'*S.U - S.V*Diagonal(S.S))/norm(Q) < 1e-9
	end

	@testset "tall-skinny and tranpose" begin
		X = randn(10000, 10)
		S1 = Severo.irlba(X, 2, tol=1e-9)
		S2 = Severo.irlba(X', 2, tol=1e-9)
		@test S1.S ≈ S2.S
		@test abs.(S1.U) ≈ abs.(S2.V) #XXX lazy testing
		@test abs.(S1.V) ≈ abs.(S2.U)
	end

	@testset "count matrix" begin
		X = sprand(3000, 10000, 0.01, (i) -> rand(Poisson(10), i))
		mu = vec(mean(X, dims=2))

		C = Severo.CenteredMatrix(X', mu)
		S = Severo.irlba(C, 10)
		SS = svd(C)
		@test S.S ≈ SS.S[1:10]
	end
end
