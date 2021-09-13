# copyright imec - evaluation license - not for distribution

using Severo, Test

@testset "metrics" begin
    @testset "purity" begin
        classes = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
        clusters = [2, 1, 2, 3, 5, 4, 5, 5, 5, 2, 5]
        @test purity(clusters, classes) ≈ 8//11
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
