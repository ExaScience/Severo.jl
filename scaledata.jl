import SparseArrays: SparseMatrixCSC, SparseVector, SparseColumnView, SparseMatrixCSCView, nonzeros, nonzeroinds, nnz, nzrange

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

function mean_var(A::SparseMatrixCSC)
  n,d = size(A)
  mu = zeros(d)
  std = zeros(d)

  for (i,a) in enumerate(eachcol(A))
    mu[i], std[i] = mean_var(a)
  end

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

function log_norm(A::SparseMatrixCSC{T}; scale_factor=1.0) where {T <: Signed}
  B = similar(A, Float64)

  for (a,b) in zip(eachcol(A), eachcol(B))
    nonzeros(b) .= log1p.(scale_factor * nonzeros(a) ./ sum(nonzeros(a)))
  end

  B
end

function filter_data(A::SparseMatrixCSC{T}; min_cells=0, min_features=0) where {T <: Signed}
  features_per_cell = vec(sum(A .> 0, dims=2))
  CI = (features_per_cell .>= min_features)

  cells_per_feature = vec(sum(A .> 0, dims=1))
  FI = (cells_per_feature .>= min_cells)

  A[CI, FI]
end
