function ann(X, k; eps=0.0)
	n,d = size(X)

	nn_index = Matrix{Int32}(undef, n, k)
	distances = Matrix{Float64}(undef, n, k)

	ccall(("FindNeighbours", "./neighbours.so"), Cvoid,
		(Ptr{Float64}, Cint, Cint, Cint, Float64, Ptr{Int32}, Ptr{Float64}),
		X, n, d, k, eps, nn_index, distances)
	nn_index, distances
end

