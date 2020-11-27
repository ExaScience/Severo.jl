function filter_features(A::CountMatrix; min_cells=0)
    cells_per_feature = vec(sum(A .> 0, dims=1))
    FI = (cells_per_feature .>= min_cells)

    A[:, FI], FI
end

function filter_cells(A::CountMatrix; min_features=0, min_feature_count=0, min_umi=0)
    features_per_cell = vec(sum(A .> min_feature_count, dims=2))
    CI = (features_per_cell .>= min_features)

    if min_umi > 0
        CI .&= (vec(sum(A, dims=2)) .> min_umi)
    end

    if !all(CI)
        A = A[CI,:]
    end

    A, CI
end

"""
    filter_features(A::NamedCountMatrix; min_cells=0)

Filter a count matrix, removing features for which the metrics fall below the given thresholds

**Arguments**:

    - `A`: the count matrix
    - `min_cells`: include features detected in at least this many cells

**Return value**:

The filtered matrix with features removed
"""
function filter_features(A::NamedCountMatrix; min_cells=0)
    counts, FI = filter_features(A.array; min_cells=min_cells)
    barcodes, features = names(A)
    NamedArray(counts, (barcodes, features[FI]), A.dimnames)
end

"""
    filter_cells(A::NamedCountMatrix; min_features=0, min_feature_count=0, min_umi=0)

Filter a labeled count matrix, removing cells for which the metrics fall below the given thresholds

**Arguments**:

    - `A`: the count matrix
    - `min_features`: include cells where at least this many features are detected
    - `min_features_count`: threshold on the count for which a feature is marked "detected"
    - `min_umi`: include cells where the total of umi counts is at least this value

**Return value**:

The filtered, labeled matrix with cells removed
"""
function filter_cells(A::NamedCountMatrix; min_features=0, min_feature_count=0, min_umi=0)
    counts, CI = filter_cells(A.array; min_features=min_features, min_feature_count=min_feature_count, min_umi=min_umi)
    barcodes, features = names(A)
    NamedArray(counts, (barcodes[CI], features), A.dimnames)
end

"""
    filter_counts(A::NamedCountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)

Filter a labeled count matrix, removing cells and features for which the metrics fall below the given threshold

First cells are removed using `filter_cells` and then features using `filter_features`. This order can be important!

**Arguments**:

    - `A`: the count matrix
    - `min_cells`: include features detected in at least this many cells
    - `min_features`: include cells where at least this many features are detected
    - `min_features_count`: threshold on the count for which a feature is marked "detected"
    - `min_umi`: include cells where the total of umi counts is at least this value

**Return value**:

The filtered, labeled matrix with cells and features removed
"""
function filter_counts(A::NamedCountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)
    counts, CI = filter_cells(A.array; min_features=min_features, min_feature_count=min_feature_count, min_umi=min_umi)
    counts, FI = filter_features(counts; min_cells=min_cells)
    barcodes, features = names(A)
    NamedArray(counts, (barcodes[CI], features[FI]), A.dimnames)
end

function percentage_features(X::NamedCountMatrix, pattern::Union{Regex, AbstractString})
    s = vec(sum(X, dims=2))
    FI = findall(x -> occursin(pattern, x), names(X,2))
    sf = vec(sum(X[:,FI], dims=2))
    sf ./ s
end
