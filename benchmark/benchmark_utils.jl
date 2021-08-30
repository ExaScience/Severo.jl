# copyright imec - evaluation license - not for distribution

using BenchmarkTools
using Severo

const utils = BenchmarkGroup()

utils["counting_sort"] = @benchmarkable Severo.counting_sort(y, 100) setup=(y = rand(1:100, 10000))
utils["make_unique"] = @benchmarkable Severo.make_unique(v) setup=(v =string.(rand('a':'e', 10000)))
