using Cell, Test
import SparseArrays: sprand, sparse
import LinearAlgebra: mul!

@testset "scaling" begin
    X = sparse(Int64[
        0  0  5  3  0  0  0  0  0  2
        0  1  0  0  0  0  0  3  0  0
        0  0  0  0  6  0  0  0  0  0
        3  0  0  0  2  0  0  0  0  0
        0  6  0  0  0  0  2  0  3  0
    ])
    C = convert_counts(X)

    S, mu = scale(C)
    @test names(S,1) == ["cell-1", "cell-2", "cell-3", "cell-4", "cell-5", "cell-6", "cell-7", "cell-8", "cell-9", "cell-10"]
    @test names(S,2) == [ "gene-1", "gene-2", "gene-3", "gene-4", "gene-5"]
    @test names(mu,1) == [ "gene-1", "gene-2", "gene-3", "gene-4", "gene-5"]

    X_mu, X_std = Cell.mean_std(X')
    @test (X_mu ./ X_std) ≈ mu
    @test Matrix((X' .- X_mu') ./ X_std') ≈ Matrix(S .- mu')
end

@testset "centered matrix dense" begin
	A = Float64[1 0 0; 0 1 0; 0 0 1; 1 0 1]
	C = Cell.CenteredMatrix(A, [0.1, 0.2, 0.3])
	@test size(C) == (4,3)
	@test eltype(C) == Float64

	Q = convert(Matrix, C)

	r = [0.991, 0.228, 0.291]
	@test C*r ≈ Q*r

	y0 = [0.820, 0.430, 0.264, 0.789]
	y = copy(y0)
	@test mul!(y, C, r, 2.0, 1.0) ≈ (2*Q*r + y0)

	@test C'*y ≈ Q'*y

	r0 = copy(r)
	@test mul!(r, C', y, 2.0, 1.0) ≈ (2*Q'*y + r0)
end

@testset "centered matrix sparse" begin
	A = sparse(Float64[1 0 0; 0 1 0; 0 0 1; 1 0 1])
	C = Cell.CenteredMatrix(A, [0.1, 0.2, 0.3])
	@test size(C) == (4,3)
	@test eltype(C) == Float64

	Q = convert(Matrix, C)

	r = [0.991, 0.228, 0.291]
	@test C*r ≈ Q*r

	y0 = [0.820, 0.430, 0.264, 0.789]
	y = copy(y0)
	@test mul!(y, C, r, 2.0, 1.0) ≈ (2*Q*r + y0)
	@test @allocated(mul!(y, C, r, 2.0, 1.0)) == 0

	@test C'*y ≈ Q'*y

	r0 = copy(r)
	@test mul!(r, C', y, 2.0, 1.0) ≈ (2*Q'*y + r0)
end
