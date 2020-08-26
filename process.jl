import SparseArrays: SparseMatrixCSC
import HDF5: h5read
import Printf: @sprintf

include("scaledata.jl")
include("irlba.jl")
include("ranksum.jl")

function read_sparse(fname, dataset="/data")
  p = h5read(fname, @sprintf("%s/p", dataset))
  i   = h5read(fname, @sprintf("%s/i", dataset))
  x = h5read(fname, @sprintf( "%s/x", dataset))
  dim = h5read(fname, @sprintf("%s/shape", dataset))
  SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
end

function read_mm10(fname, dataset="/mm10")
  p = h5read(fname, @sprintf("%s/indptr", dataset))
  i   = h5read(fname, @sprintf("%s/indices", dataset))
  x = h5read(fname, @sprintf( "%s/data", dataset))
  dim = h5read(fname, @sprintf("%s/shape", dataset))
  SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
end

function read_data(fname, dataset="/mm10")
  A = read_mm10(fname, dataset)
  lbls = h5read(fname, "/idents")
  copy(A'), lbls
end

A, lbls = read_data("/data/thaber/1M_nn.h5")
A = filter_data(A; min_cells=3, min_features=200)
B = log_norm(A, scale_factor=1e4)
println("loading done")

features = h5read("/data/thaber/1M_nn.h5", "/features/id_filtered")
C = B[:, features]
mu,std = mean_var(C)
println("irlba")
S = irlba(C, 100; center=mu, scale=std)

#@time pvals = findallmarkers(A, lbls)
