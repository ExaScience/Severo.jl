using BenchmarkTools
using Cell

import SparseArrays: sprand, sparse
import Distributions: Poisson, rand
import Random: MersenneTwister

const input = BenchmarkGroup()

X = sprand(MersenneTwister(1234), 10000, 3000, 0.01, (r,i) -> rand(r, Poisson(4), i))
input["filter_counts"] = @benchmarkable Cell.filter_counts($X; min_cells=200, min_features=3)
input["log_normalize"] = @benchmarkable Cell.log_norm($X)
input["row_normalize"] = @benchmarkable Cell.row_norm($X)
