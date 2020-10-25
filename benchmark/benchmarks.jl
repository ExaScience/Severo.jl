using BenchmarkTools
using Cell

const SUITE = BenchmarkGroup()

include("benchmark_utils.jl")
SUITE["utils"] = utils

include("benchmark_input.jl")
SUITE["input"] = input

if abspath(PROGRAM_FILE) == @__FILE__
	tune!(SUITE)
	results = run(SUITE, verbose=true, seconds=10)
	BenchmarkTools.save("benchmarks.json", results)
end
