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
import StatsBase: sample
import Statistics: mean, std

function mc_ratio(v::AbstractVector{T}, w::AbstractVector{T}, n::Int64=100) where T
	r = if n != 0 && n < length(v)
		vs = sample(v, n)
		ws = sample(w, n)
		ws ./ vs
	else
		w ./ v'
	end

	r = filter(isfinite, r)
	mu = mean(r)
	sigma = std(r, mean=mu)
	mu, sigma
end

import BenchmarkTools: ratio, Trial, TrialRatio, Parameters, params
function compare(a::Trial, b::Trial; n::Int64=100)
	ttol = max(params(a).time_tolerance, params(b).time_tolerance)
	mtol = max(params(a).memory_tolerance, params(b).memory_tolerance)
	p = Parameters(params(a); time_tolerance = ttol, memory_tolerance = mtol)
	return TrialRatio(p,
				first(mc_ratio(a.times, b.times, n)),
				first(mc_ratio(a.gctimes, b.gctimes, n)),
				ratio(a.memory, b.memory),
				ratio(a.allocs, b.allocs))
end

function compare(groups::BenchmarkGroup...; kwargs...)
	BenchmarkTools.mapvals((x...) -> compare(x...; kwargs...), groups...)
end
