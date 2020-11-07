using Cell, Test

@testset "utils" begin
    @testset "make_unique" begin
        @test Cell.make_unique(["a", "a", "a"]) == ["a", "a.1", "a.2"]
        @test Cell.make_unique(["a", "b", "a", "a.1"]) == ["a", "b", "a.2", "a.1"]
    end

    @testset "count_map" begin
        @test_throws AssertionError Cell.count_labels([0, 1, 2])

        v = [5, 4, 9, 10, 10, 5, 2, 6, 10, 7]
        @test Cell.count_labels(v) == 10
        @test Cell.count_map(v, 10) == [0, 1, 0, 1, 2, 1, 1, 0, 1, 3]
    end

    @testset "counting_sort" begin
        v = [5, 4, 9, 10, 10, 5, 2, 6, 10, 7]
        ix, counts = Cell.counting_sort(v, 10)
        @test issorted(v[ix])
        @test diff([1; counts]) == [0, 1, 0, 1, 2, 1, 1, 0, 1, 3]

        @test_throws MethodError Cell.counting_sort([1.0, 2.0], 2)
    end

    @testset "rank" begin
        @test Cell.tiedrank([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]) == [4.5, 1.5, 6.0, 1.5, 8.0, 11.0, 3.0, 10.0, 8.0, 4.5, 8.0]
        @test Cell.tiedrank([1, 1, 2, 3, 3, 4, 5, 5, 5, 6, 9]) == [1.5, 1.5, 3.0, 4.5, 4.5, 6.0, 8.0, 8.0, 8.0, 10.0, 11.0]
    end
end
