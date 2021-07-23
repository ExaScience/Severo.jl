using Severo
using Test

using Distances
import Statistics: quantile

@testset "nn" begin

@testset "euclidian" begin
    X = rand(100, 10)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(X, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, X, dims=1)
    j = map(1:size(X,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .1) > .95
end

@testset "cosine" begin
    X = rand(100, 10)

    metric = CosineDist()
    k = 4
    nn_index, distances = Severo.ann(X, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, X, dims=1)
    j = map(1:size(X,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .1) > .95
end

@testset "view" begin
    X = rand(100, 20)
    Z = view(X, :, 1:10)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(Z, k, metric)
    @test eltype(distances) == Float64

    D = pairwise(metric, Z, dims=1)
    j = map(1:size(Z,1)) do i
        nn = partialsortperm(view(D,:,i), 1:k)
        x = length(intersect(Set(nn),Set(nn_index[i,:])))
        x / (k + (k - x))
    end

    @test quantile(j, .1) > .95
end

@testset "32bit" begin
    X = rand(Float32, 100, 20)

    metric = Euclidean()
    k = 4
    nn_index, distances = Severo.ann(X, k, metric)
    @test eltype(distances) == Float32

    metric = CosineDist()
    nn_index, distances = Severo.ann(X, k, metric)
    @test eltype(distances) == Float32

    metric = Euclidean()
    Z = view(X, :, 1:10)
    nn_index, distances = Severo.ann(Z, k, metric)
    @test eltype(distances) == Float32
end

end
