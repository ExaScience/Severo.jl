# copyright imec - evaluation license - not for distribution

using BenchmarkTools
using Severo

import Random: MersenneTwister
using NamedArrays

const embedding = BenchmarkGroup()

X = NamedArray(rand(MersenneTwister(1234), 10000, 100))

embedding["pca"] = @benchmarkable Severo.embedding($X, 2, method=:pca, tol=1e-5)
embedding["umap"] = @benchmarkable Severo.embedding($X, 2, method=:umap, init=:random)
