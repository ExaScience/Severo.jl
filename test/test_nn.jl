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
using Severo
using Test

using Distances
import Statistics: quantile
import Random: default_rng

@testset "nn" begin

@testset "euclidian" begin
    X = rand(100, 10)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(default_rng(), X, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, X, dims=1)
    j = map(1:size(X,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .3) == 1.0
end

@testset "cosine" begin
    X = rand(100, 10)

    metric = CosineDist()
    k = 4
    nn_index, distances = Severo.ann(default_rng(), X, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, X, dims=1)
    j = map(1:size(X,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .3) == 1.0
end

@testset "view" begin
    X = rand(100, 20)
    Z = view(X, :, 1:10)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(default_rng(), Z, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, Z, dims=1)
    j = map(1:size(Z,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .3) == 1.0
end

@testset "32bit" begin
    X = rand(Float32, 100, 20)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(default_rng(), X, k, metric)
    @test eltype(distances) == Float32

    metric = CosineDist()
    nn_index, distances = Severo.ann(default_rng(), X, k, metric)
    @test eltype(distances) == Float32

    metric = Euclidean()
    Z = view(X, :, 1:10)
    nn_index, distances = Severo.ann(default_rng(), Z, k, metric)
    @test eltype(distances) == Float32
end

end
