using Cell, Test
import SparseArrays: sprand, sparse
import Distributions: Poisson, rand

@testset "input" begin
    X = sprand(10, 100, .1, (i) -> rand(Poisson(4), i))
    C = convert_counts(X)
    @test size(C) == (100,10)
    @test isa(C, NamedCountMatrix)

    X = sparse(Int64[
        0  0  5  3  0  0  0  0  0  2
        0  1  0  0  0  0  0  3  0  0
        0  0  0  0  6  0  0  0  0  0
        3  0  0  0  2  0  0  0  0  0
        0  6  0  0  0  0  2  0  3  0
    ])
    C = convert_counts(X)
    C = filter_counts(C; min_features=1, min_cells=2, min_umi=2)
    @test size(C) == (7, 4)
    @test names(C,1) == ["cell-1", "cell-2", "cell-3", "cell-4", "cell-5", "cell-8", "cell-9"]
    @test names(C,2) == [ "gene-1", "gene-2", "gene-4", "gene-5"]

    Q = normalize(C; method=:relativecounts)
    @test all(sum(Q, dims=2) .≈ 1.0)

    Q2 = normalize(C; method=:lognormalize)
    all(Q2 .≈ log1p.(Q))

    Q = normalize(C; method="relativecounts", scale_factor=10.)
    @test all(sum(Q, dims=2) .≈ 10.0)
end
