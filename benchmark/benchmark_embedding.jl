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

import Random: MersenneTwister
using NamedArrays

const embedding = BenchmarkGroup()

X = NamedArray(rand(MersenneTwister(1234), 10000, 100))

embedding["pca"] = @benchmarkable Severo.embedding($X, 2, method=:pca, tol=1e-5)
embedding["umap"] = @benchmarkable Severo.embedding($X, 2, method=:umap, init=:random)
