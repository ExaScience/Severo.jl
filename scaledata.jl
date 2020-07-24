import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange

import HDF5: h5read
function read_data(fname="/data/thaber/scale.h5")
	p = h5read(fname, "data/p")
	i = h5read(fname, "data/i")
	x = h5read(fname, "data/x")
	dim = h5read(fname, "data/dim")
	A = SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
	X = h5read(fname, "scaled")
	copy(A'), X
end
A, X = read_data()

function mean_var(x::Union{SparseColumnView, SparseVector})
   n = length(x)

   count = n - nnz(x)
   mu = s = zero(eltype(x))

   # nonzeros
   for v in nonzeros(x)
     count += 1
     delta = (v - mu)
     mu += delta / count
     s += delta * (v - mu)
   end

   std = sqrt(s / (n-1))
   mu, std
end

function scale_center(A::SparseMatrixCSC)
  n,d = size(A)
  B = similar(A)

  mu = zeros(d)
  for (i,(a,b)) in enumerate(zip(eachcol(A), eachcol(B)))
    mu[i], std = mean_var(a)
    nonzeros(b) .= nonzeros(a) ./ std
    mu[i] /= std
  end
  B, mu
end
