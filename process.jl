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
  features = h5read(fname, "/features/id")
  copy(A'), lbls, features
end

A, lbls, features = read_data("/data/thaber/1M_nn.h5")
A, _, FI = filter_data(A; min_cells=3, min_features=200)
B = log_norm(A, scale_factor=1e4)

@time pvals = findallmarkers(A, lbls)
