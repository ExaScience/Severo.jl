using HDF5
using LinearAlgebra
using Printf
using SparseArrays

function read_sparse(fname, dataset="/data")
	p = h5read(fname, @sprintf("%s/indptr", dataset))
	i = h5read(fname, @sprintf("%s/indices", dataset))
	x = h5read(fname, @sprintf("%s/data", dataset))
	dim = h5read(fname, @sprintf("%s/shape", dataset))
	SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
end

include("clustering.jl")

Z = h5read("/data/mca_res.h5", "/pca/embeddings")[:,1:10]
SNN = read_sparse("/data/mca_res.h5", "/snn")

assignment = modularity_cluster(SNN; resolution=0.5)
