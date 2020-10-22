using Cell, Test

@testset "utils" begin
	@testset "make_unique" begin
		@test Cell.make_unique(["a", "a", "a"]) == ["a", "a.1", "a.2"]
		@test Cell.make_unique(["a", "b", "a", "a.1"]) == ["a", "b", "a.2", "a.1"]
	end
end
