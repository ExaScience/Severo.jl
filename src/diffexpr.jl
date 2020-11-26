import DataFrames: DataFrame
import SparseArrays: SparseColumnView, SparseVector, nonzeros, nonzeroinds
import Distributions: TDist

include("ranksum.jl")

#function wilcoxon_test(x::SparseVec, lbls::AbstractVector{<:Integer}, nlabels::Int64=count_labels(lbls))
function wilcoxon_test(x::SparseVec, lbls::AbstractVector{Bool})
    ranksumtest(x, lbls)
end

function unequal_var_ttest(m1, v1, n1, m2, v2, n2)
    d = m1 - m2
    vn1 = v1 / n1
    vn2 = v2 / n2

    df = (vn1 + vn2)^2 / (vn1^2 / (n1 - 1) + vn2^2 / (n2 - 1))
    if isnan(df)
        df = 1.
    end

    t = d / sqrt(vn1 + vn2)
    prob = 2.0 * ccdf(TDist(df), abs(t))

    t, prob
end

function t_test(x::SparseVec, lbls::AbstractVector{<:Integer}, nlabels::Int64=count_labels(lbls))
    mu, var, n = mean_var(x, lbls, nlabels)
    mu2 = (sum(n .* mu) .- (n .* mu)) ./ (length(x) .- n)

    var2 = var .* (n.-1) .+ mu.^2 .*n
    var2 .= (sum(var2) .- var2 .- mu2.^2 .* (length(x) .- n)) ./ (length(x) .- n .- 1)

    scores = Vector{Float64}(undef, nlabels)
    pvals = Vector{Float64}(undef, nlabels)
    @inbounds for k in 1:nlabels
        n2 = length(x) - n[k]
        scores[k], pvals[k] = unequal_var_ttest(mu[k], var[k], n[k], mu2[k], var2[k], n2)
    end

    scores, pvals
end

function counts(ix::AbstractVector{<:Integer}, v::AbstractVector{T}, lbls::AbstractVector{<:Integer},
        nlabels::Integer=count_labels(lbls); thresh_min::T = zero(T)) where T
    n = zeros(Int64, nlabels)
    for (i,x) in zip(ix, v)
        if x > thresh_min
            n[lbls[i]] += 1
        end
    end
    n
end

function counts(x::SparseColumnView{T}, lbls::AbstractVector{<:Integer}, nlabels::Integer=length(unique(lbls)); thresh_min::T = zero(T)) where T
    counts(nonzeroinds(x), nonzeros(x), nlabels; thresh_min=thresh_min)
end

function log_means(ix::AbstractVector{<:Integer}, v::AbstractVector, lbls::AbstractVector{<:Integer},
        nlabels::Integer=count_labels(lbls), nc::AbstractVector{<:Integer}=count_map(lbls, nlabels))

    mu1 = zeros(Float64, nlabels)
    n = zeros(Int64, nlabels)
    @inbounds for (i,x) in zip(ix, v)
        idx = lbls[i]
        n[idx] += 1
        mu1[idx] += (x - mu1[idx]) / n[idx]
    end

    @inbounds for i in 1:length(n)
        mu1[i] *= n[i] / nc[i]
    end

    mu2 = (sum(nc .* mu1) .- (nc .* mu1)) ./ (length(lbls) .- nc)

    mu1 .= log1p.(mu1)
    mu2 .= log1p.(mu2)
    mu1, mu2
end

function log_means(x::SparseVec, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls),
        nc::AbstractVector{<:Integer}=count_map(lbls, nlabels))
    log_means(nonzeroinds(x), nonzeros(x), lbls, nlabels, nc)
end

function log_means_exp(x::SparseVec, lbls::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls),
        nc::AbstractVector{<:Integer}=count_map(lbls, nlabels))
    log_means(nonzeroinds(x), expm1.(nonzeros(x)), lbls, nlabels, nc)
end

"""
    filter_features(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer};
            logfc_threshold::Float64=0.0, min_pct::Float64=0.0, min_diff_pct::Float64=-Inf)

Filter features for each of the classes in a dataset.

**Arguments**:

    -``X``: count or data matrix
    -``idents``: class identity for each cell
    -``logfc_threshold``: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
    -``min_pct``: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
    -``min_diff_pct``: only test genes that show a minimum difference in the fraction of detection between the two groups.

**Return values**:

Selection matrix for each feature and class
"""
function filter_features(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer}, nlabels::Int64=count_labels(idents);
        logfc_threshold::Float64=0.25, min_pct::Float64=0.1, min_diff_pct::Float64=-Inf, log::Bool=false)
    @assert ! isa(X, NamedCountMatrix) || !log "strange combination of count data and log"

    ncells, nfeatures = size(X)
    @assert ncells == length(idents)

    # fix label order if necessary
    lbls = if names(X, 1) != names(idents, 1)
        idents[names(X,1)].array
    else
        idents.array
    end

    mean_fun = if log
        log_means_exp
    else
        log_means
    end

    # number of cells per class
    nc = count_map(lbls, nlabels)

    selected = falses(nlabels, nfeatures)
    for (i, x) in enumerate(eachcol(X.array))
        C = counts(nonzeroinds(x), nonzeros(x), lbls, nlabels)
        OC = sum(C) .- C

        pct1 = round.(C ./ nc; digits=3)
        pct2 = round.(OC ./ (ncells .- nc); digits=3)

        alpha_min = max.(pct1, pct2)
        alpha_diff = abs.(pct1 .- pct2)

        selected[:,i] = ((alpha_min .> min_pct) .& (alpha_diff .> min_diff_pct))
        any(selected[:,i]) || continue

        mu1, mu2 = mean_fun(x, lbls, nlabels, nc)
        selected[:,i] .&= abs.(mu1 .- mu2) .> logfc_threshold
    end

    label_names = map(string, 1:nlabels)
    feature_names = names(X, 2)
    NamedArray(selected, (label_names, feature_names), (:labels, dimnames(X,2)))
end

"""
    findmarkers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer}; method=:wilcoxon, kw...)

Finds markers (differentially expressed genes) for each of the classes in a dataset.

**Arguments**:

    -``X``: count or data matrix
    -``idents``: class identity for each cell
    -``method``: Which test to use, supported are: [wilcoxon]
    -``logfc_threshold``: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
    -``kw...``: additional parameters passed down to the method

**Return values**:

"""
function findmarkers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer}; method=:wilcoxon, log::Bool=false, kw...)
    @assert ! isa(X, NamedCountMatrix) || !log "strange combination of count data and log"

    if isa(method, AbstractString)
        method = Symbol(method)
    end

    f = if method == :wilcoxon
        wilcoxon_test
    elseif method == :t
        t_test
    else
        error("unknown differential expression method: $method")
    end

    ncells, nfeatures = size(X)
    @assert ncells == length(idents)

    # fix label order if necessary
    lbls = if names(X, 1) != names(idents, 1)
        idents[names(X,1)].array
    else
        idents.array
    end

    mean_fun = if log
        log_means_exp
    else
        log_means
    end

    # number of cells per class
    nlabels = count_labels(lbls)
    nc = count_map(lbls, nlabels)

    scores = Matrix{Float64}(undef, nfeatures, nlabels)
    pvals = Matrix{Float64}(undef, nfeatures, nlabels)
    logfcs = Matrix{Float64}(undef, nfeatures, nlabels)

    for (i,x) in enumerate(eachcol(X.array))
        mu1, mu2 = mean_fun(x, lbls, nlabels, nc)

        score, pval = f(x, lbls, nlabels; kw...)
        scores[i, :] = score
        pvals[i, :] = pval
        logfcs[i, :] = mu1 .- mu2
    end

    scores, pvals, logfcs
end
