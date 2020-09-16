function ann(X, k; eps=0.0)
	n,d = size(X)

	nn_index = Matrix{Int32}(undef, n, k)
	distances = Matrix{Float64}(undef, n, k)

	ccall(("FindNeighbours", "./neighbours.so"), Cvoid,
		(Ptr{Float64}, Cint, Cint, Cint, Float64, Ptr{Int32}, Ptr{Float64}),
		X, n, d, k, eps, nn_index, distances)
	nn_index, distances
end

function compute_snn(X, k; eps=0.0, prune=1/15)
	nn_index, distances = ann(X, k, eps=eps)
	nn = sparse(repeat(1:size(nnindex,1),20), vec(nnindex), ones(length(nnindex)))

	snn = nn * nn'
	f(x) = x / (k + (k - x))
	nonzeros(snn) .= f.(nonzeros(snn))
	droptol!(snn, prune)
	snn
end
