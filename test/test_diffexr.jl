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

using Test
using Severo

import SparseArrays: SparseMatrixCSC
import NamedArrays: NamedArray, names

@testset "diffex" begin
X = SparseMatrixCSC(100, 20,
    [1, 20, 40, 57, 80, 99, 109, 117, 133, 152, 162, 173, 182, 195, 210, 219, 230, 248, 269, 289, 306],
    [3, 4, 5, 10, 11, 17, 18, 26, 31, 34, 37, 38, 40, 46, 50, 54, 66, 74, 86, 2, 3, 4, 6, 8, 11, 16,
    18, 25, 34, 42, 49, 51, 65, 68, 69, 70, 74, 80, 99, 1, 9, 12, 22, 29, 33, 39, 41, 50, 66, 69, 91, 92,
    93, 95, 96, 97, 1, 2, 3, 4, 5, 6, 10, 11, 12, 24, 31, 33, 35, 39, 55, 59, 69, 70, 90, 92, 94, 96,
    98, 3, 4, 8, 10, 14, 17, 18, 19, 22, 27, 32, 47, 48, 51, 80, 82, 92, 97, 98, 27, 33, 42, 50, 52, 53,
    54, 72, 93, 97, 15, 40, 44, 45, 80, 84, 92, 97, 4, 18, 19, 36, 38, 39, 44, 46, 48, 52, 62, 67, 70, 76,
    84, 96, 1, 9, 15, 16, 20, 26, 37, 40, 53, 54, 55, 62, 69, 80, 91, 93, 94, 98, 100, 1, 6, 9, 13, 22,
    28, 40, 54, 62, 77, 2, 8, 34, 37, 40, 43, 48, 79, 81, 86, 92, 18, 21, 24, 30, 44, 45, 59, 76, 93, 7,
    15, 21, 22, 23, 41, 44, 57, 59, 84, 87, 96, 99, 5, 6, 9, 12, 15, 24, 28, 36, 39, 51, 55, 68, 75, 98,
    99, 6, 29, 39, 43, 44, 45, 69, 73, 88, 3, 4, 14, 21, 39, 40, 48, 50, 65, 68, 80, 7, 14, 17, 18, 21,
    25, 28, 29, 42, 50, 52, 78, 82, 83, 90, 94, 95, 98, 2, 3, 8, 18, 20, 23, 31, 34, 38, 41, 45, 60, 62,
    64, 71, 73, 78, 83, 91, 99, 100, 3, 6, 12, 16, 17, 18, 21, 22, 33, 37, 49, 54, 62, 65, 72, 74, 78, 84,
    92, 100, 1, 10, 11, 26, 27, 48, 54, 59, 64, 74, 78, 85, 91, 92, 94, 99, 100],
    [ 2, 2, 1, 1, 2, 2, 19, 2, 6, 13, 6, 2, 9, 4, 12, 2, 21, 10, 12, 1, 1, 1, 1, 1, 9, 4, 5, 6, 1, 7, 5, 5, 6, 2,
    1, 1, 5, 6, 4, 1, 1, 18, 2, 11, 3, 7, 2, 11, 14, 6, 4, 9, 10, 11, 4, 6, 8, 2, 3, 1, 1, 1, 1, 12, 7, 5, 3, 18,
    1, 4, 6, 7, 6, 12, 7, 1, 9, 6, 6, 1, 2, 1, 2, 1, 3, 11, 2, 5, 5, 15, 10, 19, 10, 1, 5, 8, 16, 8, 3, 4, 4, 26,
    2, 1, 6, 2, 10, 2, 8, 3, 1, 1, 11, 3, 5, 11, 13, 3, 2, 12, 5, 8, 8, 13, 4, 4, 9, 9, 3, 9, 5, 9, 3, 13, 9, 15,
       31, 11, 3, 2, 12, 14, 5, 18, 4, 2, 10, 16, 7, 8, 7, 1, 4, 6, 4, 8, 8, 3, 11, 4, 5, 1, 6, 2, 7, 1, 29, 5, 7, 7,
    3, 2, 4, 2, 9, 10, 8, 1, 5, 6, 7, 7, 3, 4, 16, 3, 6, 9, 5, 3, 3, 3, 5, 1, 5, 3, 7, 10, 7, 7, 16, 5, 13, 3,
    7, 16, 5, 17, 7, 5, 9, 3, 4, 9, 2, 5, 13, 1, 2, 16, 3, 4, 8, 4, 8, 15, 5, 2, 2, 2, 13, 6, 1, 2, 8, 3, 3, 1,
       17, 1, 5, 5, 10, 4, 10, 5, 1, 9, 1, 4, 16, 12, 17, 13, 5, 2, 2, 7, 2, 4, 2, 4, 15, 4, 5, 3, 4, 2, 7, 9, 5, 8,
    5, 7, 13, 5, 6, 5, 7, 7, 9, 9, 3, 9, 3, 1, 3, 4, 3, 9, 12, 1, 5, 9, 4, 10, 14, 6, 3, 5, 1, 5, 3, 6, 2])

lbls = Severo.rep_each([1,2], [10, 90])

X = convert_counts(X)
lbls = NamedArray(lbls, names(X,1))

sel = NamedArray(Bool[1 0 0 0 0 0 1 1 1 0 1 1 0 0 0 0 0 0 1 0; 1 0 1 1 1 0 1 0 1 1 1 1 0 0 0 1 1 0 1 0], (["1", "2"], names(X,2)))

@testset "t-test" begin
    # computed using R's t.test
    true_scores = [
        -1.5304460868698  -0.9331775615929  -2.8450731468119   0.5745765095031  -1.526391759546  -2.0760364545872  -2.3522619060766   0.1164704425837  -0.2369240462328   0.8877975346075  0.0  -2.7801417370539   0.0300908747326   0.2699556041133  -0.0204483936094   0.7680498069416  -2.3138392185418   0.0773671336178   0.3487030036407   0.2551103361614
        1.5304460868698   0.9331775615929   2.8450731468119  -0.5745765095031   1.526391759546   2.0760364545872   2.3522619060766  -0.1164704425837   0.2369240462328  -0.8877975346075  0.0   2.7801417370539  -0.0300908747326  -0.2699556041133   0.0204483936094  -0.7680498069416   2.3138392185418  -0.0773671336178  -0.3487030036407  -0.2551103361614
   ]

    de = find_markers(X, lbls, method=:t)
    @test de[!,:score] ≈ vec(true_scores) # assumes 'de' is sorted by [feature,group]

    @test de == find_markers(X, lbls, method=:t, selection=trues(size(sel)))

    de = find_markers(X, lbls, method=:t, selection=sel)
    vs = vec(sel)
    @test de[!,:score] ≈ vec(true_scores)[vs] # assumes 'de' is sorted by [feature,group]
    @test de[!,:feature] == repeat(names(X,2), inner=2)[vs] # assumes 'de' is sorted by [feature,group]
    @test de[!,:group] == repeat([1,2], outer=size(true_scores,2))[vs] # assumes 'de' is sorted by [feature,group]
end

@testset "wilcoxon" begin
    # computed using R's wilcox.test
    true_pvals = [0.18191295173714, 0.04922813607314, 1.00000000000000, 0.00192017900493,
        0.17127188722667, 0.27459335840573, 0.33467860916225, 0.68539804935004,
        0.97322000673841, 0.03601086837360, 0.38577520003606, 0.30304142891167,
        0.83641849963529, 0.23642583425188, 0.90786998988961, 0.33580029036544,
        0.45551138123867, 0.53433311634786, 0.88879550022386, 0.79225415168002
    ]

    de = Severo.find_markers(X, lbls, method=:wilcoxon)
    @test de[!,:pval] ≈ repeat(true_pvals, inner=2) # assumes 'de' is sorted by [feature,group]

    @test de == find_markers(X, lbls, method=:wilcoxon, selection=trues(size(sel)))

    de = find_markers(X, lbls, method=:wilcoxon, selection=sel)
    vs = vec(sel)
    @test de[!,:pval] ≈ repeat(true_pvals, inner=2)[vs] # assumes 'de' is sorted by [feature,group]
    @test de[!,:feature] == repeat(names(X,2), inner=2)[vs] # assumes 'de' is sorted by [feature,group]
    @test de[!,:group] == repeat([1,2], outer=length(true_pvals))[vs] # assumes 'de' is sorted by [feature,group]
end
end
