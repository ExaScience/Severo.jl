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
import SparseArrays: sprand, sparse
import Distributions: Poisson, rand

@testset "filter" begin
    X = sprand(100, 10, .1, (i) -> rand(Poisson(4), i))
    C = convert_counts(X)
    @test size(C) == (100,10)
    @test isa(C, NamedCountMatrix)

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
    C = filter_counts(C; min_features=1, min_cells=2, min_umi=2)

    @test size(C) == (7, 4)
    @test names(C,1) == ["cell-1", "cell-2", "cell-3", "cell-4", "cell-5", "cell-8", "cell-9"]
    @test names(C,2) == [ "gene-1", "gene-2", "gene-4", "gene-5"]
end

@testset "normalize_cells" begin
    @testset "normalize64" begin
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
        C = filter_counts(C; min_features=1, min_cells=2, min_umi=2)

        Q = normalize_cells(C; method=:relativecounts)
        @test eltype(Q) == Float64
        @test all(sum(Q, dims=2) .≈ 1.0)
        @test names(Q,1) == ["cell-1", "cell-2", "cell-3", "cell-4", "cell-5", "cell-8", "cell-9"]
        @test names(Q,2) == [ "gene-1", "gene-2", "gene-4", "gene-5"]

        Q2 = normalize_cells(C; method=:lognormalize)
        @test eltype(Q2) == Float64
        all(Q2 .≈ log1p.(Q))

        Q = normalize_cells(C; method="relativecounts", scale_factor=10.)
        @test eltype(Q) == Float64
        @test all(sum(Q, dims=2) .≈ 10.0)
    end

    @testset "normalize32" begin
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
        C = filter_counts(C; min_features=1, min_cells=2, min_umi=2)

        Q = normalize_cells(C; method=:relativecounts, dtype=Float32)
        @test eltype(Q) == Float32

        Q2 = normalize_cells(C; method=:lognormalize, dtype=Float32)
        @test eltype(Q2) == Float32

        Q = normalize_cells(C; method="relativecounts", scale_factor=10., dtype=Float32)
        @test eltype(Q) == Float32
    end
end
