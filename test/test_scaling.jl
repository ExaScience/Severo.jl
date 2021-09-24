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
import SparseArrays: sprand, sparse, SparseVector
import LinearAlgebra: mul!
import Statistics

@testset "scaling" begin
    @testset "scaling64" begin
        X = sparse(Int64[
            0  0  0  3  0
            0  1  0  0  6
            5  0  0  0  0
            3  0  0  0  0
            0  0  6  2  0
            0  0  0  0  0
            0  0  0  0  2
            0  3  0  0  0
            0  0  0  0  3
            2  0  0  0  0
        ])
        C = convert_counts(X)

        S = scale_features(C)
        @test eltype(S) == Float64
        @test names(S,1) == ["cell-1", "cell-2", "cell-3", "cell-4", "cell-5", "cell-6", "cell-7", "cell-8", "cell-9", "cell-10"]
        @test names(S,2) == [ "gene-1", "gene-2", "gene-3", "gene-4", "gene-5"]

        X_mu, X_std = Severo.mean_std(X)
        @test eltype(X_mu) == Float64 && eltype(X_std) == Float64
        @test (X_mu ./ X_std) ≈ S.mu
        @test Matrix((X .- X_mu') ./ X_std') ≈ convert(Matrix, S)
    end

    @testset "scaling32" begin
        X = sparse(Int64[
            0  0  0  3  0
            0  1  0  0  6
            5  0  0  0  0
            3  0  0  0  0
            0  0  6  2  0
            0  0  0  0  0
            0  0  0  0  2
            0  3  0  0  0
            0  0  0  0  3
            2  0  0  0  0
        ])
        C = convert_counts(X)

        S = scale_features(C, dtype=Float32)
        @test eltype(S) == Float32

        X_mu, X_std = Severo.mean_std(Float32, X)
        @test eltype(X_mu) == Float32 && eltype(X_std) == Float32
    end
end

@testset "centered matrix" begin
    @testset "dense" begin
        A = Float64[1 0 0; 0 1 0; 0 0 1; 1 0 1]
        C = Severo.CenteredMatrix(A, [0.1, 0.2, 0.3])
        @test size(C) == (4,3)
        @test eltype(C) == Float64

        Q = convert(Matrix, C)

        r = [0.991, 0.228, 0.291]
        @test C*r ≈ Q*r

        y0 = [0.820, 0.430, 0.264, 0.789]
        y = copy(y0)
        @test mul!(y, C, r, 2.0, 1.0) ≈ (2*Q*r + y0)

        @test C'*y ≈ Q'*y

        r0 = copy(r)
        @test mul!(r, C', y, 2.0, 1.0) ≈ (2*Q'*y + r0)
    end

    @testset "sparse" begin
        A = sparse(Float64[1 0 0; 0 1 0; 0 0 1; 1 0 1])
        C = Severo.CenteredMatrix(A, [0.1, 0.2, 0.3])
        @test size(C) == (4,3)
        @test eltype(C) == Float64

        Q = convert(Matrix, C)

        r = [0.991, 0.228, 0.291]
        @test C*r ≈ Q*r

        y0 = [0.820, 0.430, 0.264, 0.789]
        y = copy(y0)
        @test mul!(y, C, r, 2.0, 1.0) ≈ (2*Q*r + y0)
        @test @allocated(mul!(y, C, r, 2.0, 1.0)) == 0

        @test C'*y ≈ Q'*y

        r0 = copy(r)
        @test mul!(r, C', y, 2.0, 1.0) ≈ (2*Q'*y + r0)
    end
end

@testset "mean_var sparse" begin
    x = SparseVector(100,
        [3, 4, 5, 10, 11, 17, 18, 26, 31, 34, 37, 38, 40, 46, 50, 54, 66, 74, 86],
        [2, 2, 1, 1, 2, 2, 19, 2, 6, 13, 6, 2, 9, 4, 12, 2, 21, 10, 12])

    lbls, nlabels = Severo.rep_each([1,2], [10, 90]), 2
    mu, var, n = Severo.mean_var(x, lbls, nlabels)
    @test n ≈ [count(lbls.==i) for i in 1:nlabels]
    @test mu ≈ [Statistics.mean(x[lbls.==i]) for i in 1:nlabels]
    @test var ≈ [Statistics.var(x[lbls.==i]) for i in 1:nlabels]

    lbls, nlabels = rand(1:4, 100), 4
    mu, var, n = Severo.mean_var(x, lbls, nlabels)
    @test n ≈ [count(lbls.==i) for i in 1:nlabels]
    @test mu ≈ [Statistics.mean(x[lbls.==i]) for i in 1:nlabels]
    @test var ≈ [Statistics.var(x[lbls.==i]) for i in 1:nlabels]
end
