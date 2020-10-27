using BenchmarkTools
using Cell

import Arpack: svds
import SparseArrays: sprand, sparse
import Distributions: Poisson, rand, randn
import Random: MersenneTwister

const irlba = BenchmarkGroup()

X = sprand(MersenneTwister(1234), 10000, 3000, 0.01, (r,i) -> rand(r, Poisson(4), i))
mu = vec(mean(X, dims=2))
C = Cell.CenteredMatrix(X', mu)

X = randn(1000, 500)

irlba["irlba_dense"] = @benchmarkable Cell.irlba($X, 50, tol=1e-5)
irlba["irlba_sparse"] = @benchmarkable Cell.irlba($C, 10, tol=1e-5)
irlba["arpack_dense"] = @benchmarkable svds($X; nsv=50, tol=1e-5)
irlba["arpack_sparse"] = @benchmarkable svds($C; nsv=10, tol=1e-5)
