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
