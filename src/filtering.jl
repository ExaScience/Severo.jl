"""
    filter_counts(A::CountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)

Filter a count matrix, removing cells and features for which the metrics fall below the given threshold

**Arguments**:

    - `A`: the count matrix
    - `min_cells`: include features detected in at least this many cells
    - `min_features`: include cells where at least this many features are detected
    - `min_features_count`: threshold on the count for which a feature is marked "detected"
    - `min_umi`: include cells where the total of umi counts is at least this value

**Return value**:

The filtered matrix with cells and features removed
"""
function filter_counts(A::CountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)
    features_per_cell = vec(sum(A .> min_feature_count, dims=2))
    CI = (features_per_cell .>= min_features)

    if min_umi > 0
        CI .&= (vec(sum(A, dims=2)) .> min_umi)
    end

    if !all(CI)
        A = A[CI,:]
    end

    cells_per_feature = vec(sum(A .> 0, dims=1))
    FI = (cells_per_feature .>= min_cells)

    A[:, FI], CI, FI
    #  A[CI, FI] # faster but slightly different
end

"""
    filter_counts(A::NamedCountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)

Filter a labeled count matrix, removing cells and features for which the metrics fall below the given threshold

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
    counts, CI, FI = filter_counts(A.array; min_cells=min_cells, min_features=min_features, min_feature_count=min_feature_count, min_umi=min_umi)
    barcodes, features = names(A)
    NamedArray(counts, (barcodes[CI], features[FI]), (:cells, :features))
end

function percentage_features(X::NamedCountMatrix, pattern::Union{Regex, AbstractString})
    s = vec(sum(X, dims=2))
    FI = findall(x -> occursin(pattern, x), names(X,2))
    sf = vec(sum(X[:,FI], dims=2))
    sf ./ s
end
