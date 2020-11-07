using BenchmarkTools
using Cell

const SUITE = BenchmarkGroup()

include("benchmark_utils.jl")
SUITE["utils"] = utils

include("benchmark_input.jl")
SUITE["input"] = input

include("benchmark_irlba.jl")
SUITE["irlba"] = irlba

include("benchmark_embedding.jl")
SUITE["embedding"] = embedding

if abspath(PROGRAM_FILE) == @__FILE__
	tune!(SUITE)
	results = run(SUITE, verbose=true, seconds=10)

	if isfile("benchmarks.json")
		prev = BenchmarkTools.load("benchmarks.json")[1]
		include("utils.jl")
		c = compare(prev, results, n=1000)
		show(judge(c))
	end

	BenchmarkTools.save("benchmarks.json", results)
end
