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

using Severo, Test

@testset verbose=true "all" begin
include("test_utils.jl")
include("test_input.jl")
include("test_scaling.jl")
include("test_cluster.jl")
include("test_ranksum.jl")
include("test_nn.jl")
include("test_irlba.jl")
include("test_metrics.jl")
end
