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

const utils = BenchmarkGroup()

utils["counting_sort"] = @benchmarkable Severo.counting_sort(y, 100) setup=(y = rand(1:100, 10000))
utils["make_unique"] = @benchmarkable Severo.make_unique(v) setup=(v =string.(rand('a':'e', 10000)))
