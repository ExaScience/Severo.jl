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
using LinearAlgebra

@testset "metrics" begin
    @testset "purity" begin
        classes = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
        clusters = [2, 1, 2, 3, 5, 4, 5, 5, 5, 2, 5]
        @test purity(clusters, classes) ≈ 8//11
    end

    @testset "agreement" begin
        X = randn(2000, 20)
        # agreement(X, X) should be high

        samples = map(1:20) do i
            Y = randn(2000, 20)
            agreement(X, X, Y) # agreement should be low between X and Y
        end
        q = quantile(map((x)-> -(x...), samples), (0.025, 0.975))
        @test !(first(q) < 0.0 < last(q))

        samples = map(1:20) do i
            # random rotation
            Q, _ = qr(randn(20))
            Z = X * Q
            agreement(X, X, Z) # should be similar
        end
        q = quantile(map((x)-> -(x...), samples), (0.025, 0.975))
        @test (first(q) < 0.0 < last(q))
    end

    @testset "alignment" begin
        datasets = [1:1000, 1001:2000]

        X = randn(2000, 20)
        align = alignment(X, datasets...)
        # alignment should be high
        @test mean(align) > .95

        X = vcat(randn(1000, 20) .+ 10, randn(1000, 20) .- 10)
        align = alignment(X, datasets...)
        # alignment should be low
        @test mean(align) < 0.01
    end
end
