# copyright imec - evaluation license - not for distribution

using BenchmarkTools
using Cell

const utils = BenchmarkGroup()

utils["counting_sort"] = @benchmarkable Cell.counting_sort(y, 100) setup=(y = rand(1:100, 10000))
utils["make_unique"] = @benchmarkable Cell.make_unique(v) setup=(v =string.(rand('a':'e', 10000)))
