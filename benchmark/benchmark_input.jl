# copyright imec - evaluation license - not for distribution

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
