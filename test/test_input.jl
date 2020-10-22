using Cell, Test
import SparseArrays: sprand
import Distributions: Poisson, rand

@testset "input" begin
	X = sprand(10, 100, .1, (i) -> rand(Poisson(4), i))
	C = convert_data(X)
	@test size(C) == (100,10)
end
