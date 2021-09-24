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

import SparseArrays: sprand, sparse
import Distributions: Poisson, rand
import Random: MersenneTwister

const input = BenchmarkGroup()

X = sprand(MersenneTwister(1234), 10000, 3000, 0.01, (r,i) -> rand(r, Poisson(4), i))
input["filter_counts"] = @benchmarkable Severo.filter_counts($X; min_cells=200, min_features=3)
input["log_normalize"] = @benchmarkable Severo.log_norm($X)
input["row_normalize"] = @benchmarkable Severo.row_norm($X)
