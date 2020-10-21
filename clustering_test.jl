using HDF5
using LinearAlgebra
using Printf
using SparseArrays
using Random

function read_sparse(fname, dataset="/data")
	p = h5read(fname, @sprintf("%s/indptr", dataset))
	i = h5read(fname, @sprintf("%s/indices", dataset))
	x = h5read(fname, @sprintf("%s/data", dataset))
	dim = h5read(fname, @sprintf("%s/shape", dataset))
	SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
end

include("clustering.jl")
include("neighbours.jl")
include("modularity.jl")

#Z = h5read("/data/mca_res.h5", "/pca/embeddings")[:,1:10]
#snn = compute_snn(Z, 20)
snn = read_sparse("/data/mca_res.h5", "/snn")

@time assignment = modularity_cluster(snn; resolution=1.0, nrandomstarts=1)
n = Network(snn)
Random.seed!(123456)
@time cl = louvain(n; resolution=1.0, min_modularity=0.0)
@assert modularity(cl) ≈ 0.8943628471615616
