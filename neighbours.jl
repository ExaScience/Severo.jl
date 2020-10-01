import SparseArrays: sparse, nonzeros, droptol!

function ann(X, k; ntables=2*size(X,2))
	n,d = size(X)

	nn_index = Matrix{Int32}(undef, n, k)
	distances = Matrix{Float64}(undef, n, k)

	ccall(("FindNeighbours", "./neighbours.so"), Cvoid,
		(Ptr{Float64}, Cint, Cint, Cint, Cint, Ptr{Int32}, Ptr{Float64}),
		X, n, d, k, ntables, nn_index, distances)
	nn_index, distances
end

function compute_snn(X, k; ntables=2*size(X,2), prune=1/15)
	nn_index, distances = ann(X, k, ntables=ntables)
	nn = sparse(repeat(1:size(nn_index,1),k), vec(nn_index), ones(length(nn_index)))

	snn = nn * nn'
	f(x) = x / (k + (k - x))
	nonzeros(snn) .= f.(nonzeros(snn))
	droptol!(snn, prune)
	snn
end
