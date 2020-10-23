using BenchmarkTools
using Cell

const SUITE = BenchmarkGroup()

include("benchmark_utils.jl")
SUITE["utils"] = utils

#=
tune!(SUITE)
results = run(SUITE, verbose=true, seconds=10)
BenchmarkTools.save("benchmarks.json", results)
=#
