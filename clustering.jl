import SparseArrays: SparseMatrixCSC, nonzeros, getcolptr, rowvals, nnz

function modularity_cluster(SNN::SparseMatrixCSC; modularity=1, resolution=0.8, algorithm=1, randomseed=0, nrandomstarts=10, niterations=10, verbose=true)
	m,n = size(SNN)

	assignment = Vector{Int32}(undef, max(m,n));
	ccall(("ModularityClustering", "./clustering.so"), Cvoid,
		(Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Int32, Int32, Int64,
		Int32, Float64, Int32, Int32, Int32, Int32, Bool, Ptr{Int32}),
		getcolptr(SNN), rowvals(SNN), nonzeros(SNN), m, n, nnz(SNN),
		modularity, resolution, algorithm, nrandomstarts, niterations,
		randomseed, verbose, assignment)

	assignment .+= 1
end

