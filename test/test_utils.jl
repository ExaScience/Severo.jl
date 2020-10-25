using Cell, Test

@testset "utils" begin
    @testset "make_unique" begin
        @test Cell.make_unique(["a", "a", "a"]) == ["a", "a.1", "a.2"]
        @test Cell.make_unique(["a", "b", "a", "a.1"]) == ["a", "b", "a.2", "a.1"]
    end

    @testset "counting_sort" begin
        v = [5, 4, 9, 10, 10, 5, 2, 6, 10, 7]
        ix, counts = Cell.counting_sort(v, 10)
        @test issorted(v[ix])
        @test diff([1; counts]) == [0, 1, 0, 1, 2, 1, 1, 0, 1, 3]
    end
end
