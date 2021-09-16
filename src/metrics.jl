# copyright imec - evaluation license - not for distribution
import StatsBase: counts, IntegerArray, sample
import Random: AbstractRNG, default_rng

"""
    purity(clusters::IntegerArray, classes::IntegerArray)

Calculates purity between clusters and external clustering (true clusters/classes).

**Arguments**:

    - `clusters`: clustering for which to calculate purity
    - `classes`: clustering/classes with which to compare

**Return values**:

Purity score in the range [0, 1], with a score of 1 representing a pure/accurate clustering
"""
function purity(clusters::IntegerArray, classes::IntegerArray)
    c = counts(classes, clusters)
    sum(maximum(c, dims=1)) / length(classes)
end

_jaccard_index(a::IntegerArray, b::IntegerArray) = (i = length(intersect(a,b)); i / (length(a) + length(b) - i))

function _agreement(rng::AbstractRNG, nn_X::IntegerArray, Y::AbstractMatrix, k::Int64, metric::SemiMetric)
    nn_Y, _ = ann(rng, Y, k, metric)
    mean(map(_jaccard_index, eachrow(nn_X), eachrow(nn_Y)))
end

_agreement(rng::AbstractRNG, nn_X::IntegerArray, Y::NamedArray, k::Int64, metric::SemiMetric) = _agreement(rng, nn_X, Y.array, k, metric)
_agreement(rng::AbstractRNG, nn_X::IntegerArray, Y::LinearEmbedding, k::Int64, metric::SemiMetric) = _agreement(rng, nn_X, Y.coordinates, k, metric)

"""
    agreement(rng::AbstractRNG, X::AbstractMatrix, Y::Union{AbstractMatrix, LinearEmbedding}; k::Int64)

A metric for quantifying how much a transformation/factorization distorts the geometry of the
original dataset. The greater the agreement, the less distortion of geometry there is.

This is calculated by performing dimensionality reduction on the original and transformed dataset,
and measuring similarity between the k nearest neighbors for each cell in the datasets.
The Jaccard index is used to quantify similarity, and is the final metric averages across all cells.

**Arguments**:

    - `rng`: random number generator used for k-NN
    - `X`: low dimensional embedding for reference dataset
    - `Y` low dimensional embedding for transformed dataset
    - `k`: number of neighbours to find (default=15)

**Return values**:
The agreement score
"""
function agreement(rng::AbstractRNG, X::AbstractMatrix, Y::Union{AbstractMatrix, LinearEmbedding}; k::Int64=15, metric::SemiMetric=Euclidean())
    nn_X, _ = ann(rng, X, k, metric)
    _agreement(rng, nn_X, Y, k, metric)
end

"""
    agreement(rng::AbstractRNG, X::Union{AbstractMatrix, LinearEmbedding}, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64)

A metric for quantifying how much a transformation/factorization distorts the geometry of the
original dataset. The greater the agreement, the less distortion of geometry there is.

This is calculated by performing dimensionality reduction on the original and transformed dataset,
and measuring similarity between the k nearest neighbors for each cell in the datasets.
The Jaccard index is used to quantify similarity, and is the final metric averages across all cells.

**Arguments**:

    - `rng`: random number generator used for k-NN
    - `X`: low dimensional embedding for reference dataset
    - `Ys` low dimensional embedding for transformed datasets
    - `k`: number of neighbours to find (default=15)

**Return values**:
The agreement score
"""
function agreement(rng::AbstractRNG, X::AbstractMatrix, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64=15, metric::SemiMetric=Euclidean())
    nn_X, _ = ann(rng, X, k, metric)
    map(Ys) do Y
        _agreement(rng, nn_X, Y, k, metric)
    end
end

function agreement(rng::AbstractRNG, X::NamedArray, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64=15, metric::SemiMetric=Euclidean())
    agreement(rng, X.array, Ys...; k=k, metric=metric)
end

function agreement(rng::AbstractRNG, X::LinearEmbedding, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64=15, metric::SemiMetric=Euclidean())
    agreement(rng, X.coordinates, Ys...; k=k, metric=metric)
end

"""
    agreement(X::Union{AbstractMatrix, LinearEmbedding}, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64)

A metric for quantifying how much a transformation/factorization distorts the geometry of the
original dataset. The greater the agreement, the less distortion of geometry there is.

This is calculated by performing dimensionality reduction on the original and transformed dataset,
and measuring similarity between the k nearest neighbors for each cell in the datasets.
The Jaccard index is used to quantify similarity, and is the final metric averages across all cells.

**Arguments**:

    - `X`: low dimensional embedding for reference dataset
    - `Y...`: low dimensional embedding for transformed dataset(s)
    - `k`: number of neighbours to find (default=15)

**Return values**:
The agreement score
"""
agreement(X::Union{AbstractMatrix, LinearEmbedding}, Ys::Union{AbstractMatrix, LinearEmbedding}...; k::Int64=15, metric::SemiMetric=Euclidean()) = agreement(default_rng(), X, Ys...; k=k, metric=metric)

function subsample_datasets(rng::AbstractRNG, min_cells, datasets)
    samples = similar(first(datasets), min_cells * length(datasets))

    for (i, ds) in enumerate(datasets)
        idx = 1 + (i - 1) * min_cells
        samples[idx:idx+min_cells-1] = sample(rng, ds, min_cells; replace=false)
    end

    samples
end

"""
    alignment(rng::AbstractRNG, X::AbstractMatrix, datasets::AbstractVector{T}...; k::Union{Nothing,Int64}=nothing) where T

Calculates the `alignment score` as defined by Butler 2018 [doi: 10.1038/nbt.4096].
It's a quantitative metric for the alignment of datasets and calculated as follows:

    1. Randomly downsample the datasets to have the same number of cells as the smallest dataset
    2. Construct a nearest-neighbor graph based on the cells’ embedding in some low dimensional space `X`.
    3. For every cell, calculate how many of its k nearest-neighbors belong to the same dataset and average this over all cells.
    4. We then normalize by the expected number of same dataset cells and scale to range from 0 to 1.

If the datasets are well-aligned, we would expect that each cells’ nearest neighbors would be evenly shared across all datasets.

**Arguments**:

    - `rng`: random number generator used by downsampling and k-NN
    - `X`: low dimensional embedding
    - `datasets`: the split into datasets
    - `k`: number of neighbours to find. By default: 1% of the total number of cells, capped by a minimum of 10 and total number of samples drawn

**Return values**:
The alignment score
"""
function alignment(rng::AbstractRNG, X::AbstractMatrix, datasets::AbstractVector{T}...; k::Union{Nothing,Int64}=nothing) where T
    N = length(datasets)
    if N == 1
        @warn "calculating alignment for 1 dataset"
        return 1.0
    end

    min_cells = minimum(length, datasets)

    # downsample each dataset to have `min_cells` cells
    sampled_cells = subsample_datasets(rng, min_cells, datasets)

    if k === nothing
        num_cells = size(X,1)
        num_sampled = length(sampled_cells)
        k = clamp(floor(Int, 0.01*num_cells), 10, num_sampled - 1)
    end

    # calculate k-NN for each cell
    Xs = X[sampled_cells,:]
    nn_index, _ = ann(rng, Xs, k, include_self=false)

    map(1:N) do i
        idx = 1 + (i - 1) * min_cells
        r = idx:idx+min_cells-1
        x̄ = mean(count(in(r), view(nn_index, r, :); dims=2))

        1 - (x̄ - k/N) / (k - k/N)
    end
end

function alignment(rng::AbstractRNG, X::NamedArray, datasets::AbstractVector{T}...; k::Union{Nothing,Int64}=nothing) where T
    alignment(rng, X.array, datasets...; k=k)
end

"""
    alignment(X::AbstractMatrix, datasets::AbstractVector{T}...; k::Union{Nothing,Int64}=nothing) where T

Calculates the `alignment score` as defined by Butler 2018 [doi: 10.1038/nbt.4096].
It's a quantitative metric for the alignment of datasets and calculated as follows:

    1. Randomly downsample the datasets to have the same number of cells as the smallest dataset
    2. Construct a nearest-neighbor graph based on the cells’ embedding in some low dimensional space `X`.
    3. For every cell, calculate how many of its k nearest-neighbors belong to the same dataset and average this over all cells.
    4. We then normalize by the expected number of same dataset cells and scale to range from 0 to 1.

If the datasets are well-aligned, we would expect that each cells’ nearest neighbors would be evenly shared across all datasets.

**Arguments**:

    - `X`: low dimensional embedding of the aligned datasets
    - `datasets`: the split into datasets
    - `k`: number of neighbours to find. By default: 1% of the total number of cells, capped by a minimum of 10 and total number of samples drawn

**Return values**:
The alignment score
"""
alignment(X::AbstractMatrix, datasets::AbstractVector{T}...; k::Union{Nothing,Int64}=nothing) where T = alignment(default_rng(), X, datasets...; k=k)
