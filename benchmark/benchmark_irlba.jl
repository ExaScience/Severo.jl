# Severo: a software package for analysis and exploration of single-cell RNA-seq datasets.
# Copyright (c) 2021 imec vzw.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, and Additional Terms
# (see below).

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.

using BenchmarkTools
using Severo

import Arpack: svds
import SparseArrays: sprand, sparse
import Distributions: Poisson, rand, randn
import Random: MersenneTwister

const irlba = BenchmarkGroup()

X = sprand(MersenneTwister(1234), 10000, 3000, 0.01, (r,i) -> rand(r, Poisson(4), i))
mu = vec(mean(X, dims=1))
C = Severo.CenteredMatrix(X, mu)

X = randn(1000, 500)

irlba["irlba_dense"] = @benchmarkable Severo.irlba($X, 50, tol=1e-5)
irlba["irlba_sparse"] = @benchmarkable Severo.irlba($C, 10, tol=1e-5)
irlba["arpack_dense"] = @benchmarkable svds($X; nsv=50, tol=1e-5)
irlba["arpack_sparse"] = @benchmarkable svds($C; nsv=10, tol=1e-5)
