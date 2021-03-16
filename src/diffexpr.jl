# copyright imec - evaluation license - not for distribution

import DataFrames: DataFrame, sort!, select!, groupby, combine, Not
import SparseArrays: nonzeros, nonzeroinds
import Distributions: TDist

include("ranksum.jl")

function wilcoxon_test!(scores::AbstractVector{Float64}, pvals::AbstractVector{Float64},
        x::SparseVec, lbls::AbstractVector{<:Integer}, sel::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls))
    s = BitArray(undef, length(lbls))

    idx = 1
    @inbounds for k in sel
        s .= lbls .== k
        scores[idx], pvals[idx] = ranksumtest(x, s)
        idx += 1
    end

    scores, pvals
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

function t_test!(scores::AbstractVector{Float64}, pvals::AbstractVector{Float64},
        x::SparseVec, lbls::AbstractVector{<:Integer}, sel::AbstractVector{<:Integer}, nlabels::Integer=count_labels(lbls))
    mu, var, n = mean_var(x, lbls, nlabels)
    mu2 = (sum(n .* mu) .- (n .* mu)) ./ (length(x) .- n)

    var2 = var .* (n.-1) .+ mu.^2 .*n
    var2 .= (sum(var2) .- var2 .- mu2.^2 .* (length(x) .- n)) ./ (length(x) .- n .- 1)

    idx = 1
    @inbounds for k in sel
        n2 = length(x) - n[k]
        scores[idx], pvals[idx] = unequal_var_ttest(mu[k], var[k], n[k], mu2[k], var2[k], n2)
        idx += 1
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
        mu1[i] *= if nc[i] != 0
            n[i] / nc[i]
        else
            0.0
        end
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
    prefilter_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer};
            logfc_threshold::Real=0.0, min_pct::Real=0.0, min_diff_pct::Real=-Inf, only_pos:Bool=false, log::Bool=false)

Filter features for each of the classes in a dataset.

**Arguments**:

    -`X`: count or data matrix
    -`idents`: class identity for each cell
    -`logfc_threshold`: Limit testing to features which show, on average, at least X-fold difference (log-scale) between the two groups of cells
    -`min_pct`: only test features that are detected in a minimum fraction of `min_pct` cells in either of the two populations
    -`min_diff_pct`: only test features that show a minimum difference in the fraction of detection between the two groups.
    -`only_pos`: only return features with positive log fold-change
    -`log`: the data is in log-scale (default = false)

**Return values**:

Selection matrix for each feature and class
"""
function prefilter_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer}, nlabels::Int64=count_labels(idents);
        logfc_threshold::Real=0.25, min_pct::Real=0.1, min_diff_pct::Real=-Inf, only_pos::Bool=false, log::Bool=false)
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

        selected[:,i] .&= if only_pos
            (mu1 .- mu2) .> logfc_threshold
        else
            abs.(mu1 .- mu2) .> logfc_threshold
        end
    end

    label_names = map(string, 1:nlabels)
    feature_names = names(X, 2)
    NamedArray(selected, (label_names, feature_names), (:labels, dimnames(X,2)))
end

"""
    find_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer};
        method=:wilcoxon, selection::Union{Nothing, NamedArray{Bool, 2}, AbstractArray{Bool,2}}=nothing, log::Bool=false, kw...)

Finds markers (differentially expressed genes) for each of the classes in a dataset.

**Arguments**:

    -`X`: count or data matrix
    -`idents`: class identity for each cell
    -`method`: Which test to use, supported are: [wilcoxon, t]
    -`selection`: a selection of features and groups that should be considered
    -`log`: the data is in log-scale (default = false)
    -`kw...`: additional parameters passed down to the method

**Return values**:

A `DataFrame` containing a list of putative markers with associated statistics (p-values and scores) and log fold-changes.
"""
function find_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer};
        method=:wilcoxon, log::Bool=false, selection::Union{Nothing, NamedArray{Bool, 2}, AbstractArray{Bool,2}}=nothing, kw...)
    @assert ! isa(X, NamedCountMatrix) || !log "strange combination of count data and log"

    if isa(method, AbstractString)
        method = Symbol(method)
    end

    f! = if method == :wilcoxon
        wilcoxon_test!
    elseif method == :t
        t_test!
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

    # fix label order if necessary
    sel = if selection isa NamedArray{Bool,2}
         convert(BitArray, if names(X, 2) != names(selection, 2)
                selection[names(X,2)].array
            else
                selection.array
            end)
    else
        selection
    end

    mean_fun = if log
        log_means_exp
    else
        log_means
    end

    find_markers(X, lbls, f!, mean_fun, sel, kw)
end

function find_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, lbls::AbstractVector{<:Integer}, f!::Function, mean_fun::Function, ::Nothing, kw)
    ncells, nfeatures = size(X)

    # number of cells per class
    nlabels = count_labels(lbls)
    nc = count_map(lbls, nlabels)

    scores = Matrix{Float64}(undef, nlabels, nfeatures)
    pvals = Matrix{Float64}(undef, nlabels, nfeatures)
    logfc = Matrix{Float64}(undef, nlabels, nfeatures)

    @inbounds for (i,x) in enumerate(eachcol(X.array))
        mu1, mu2 = mean_fun(x, lbls, nlabels, nc)

        f!(view(scores, :,i), view(pvals, :,i), x, lbls, 1:nlabels, nlabels; kw...)
        logfc[:, i] = mu1 .- mu2
    end

    DataFrame(score=vec(scores), pval=vec(pvals), logfc=vec(logfc),
              group=repeat(1:nlabels, outer=nfeatures), feature=repeat(names(X,2), inner=nlabels))
end

function find_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, lbls::AbstractVector{<:Integer}, f!::Function,
        mean_fun::Function, selection::BitArray{2}, kw)
    ncells, nfeatures = size(X)

    # number of cells per class
    nlabels = count_labels(lbls)
    nc = count_map(lbls, nlabels)

    nf = vec(count(selection, dims=1))
    n = sum(nf)

    scores = Vector{Float64}(undef, n)
    pvals = Vector{Float64}(undef, n)
    logfc = Vector{Float64}(undef, n)
    groups = Vector{Int64}(undef, n)
    features = Vector{String}(undef, n)

    feature_names = names(X,2)

    idx = 1
    @inbounds for i in 1:nfeatures
        nf[i] == 0 && continue

        x = view(X.array, :, i)
        s = findall(view(selection, :, i))
        r = idx:idx+nf[i]-1

        mu1, mu2 = mean_fun(x, lbls, nlabels, nc)
        logfc[r] .= (mu1 .- mu2)[s]

        f!(view(scores, r), view(pvals, r), x, lbls, s, nlabels; kw...)

        groups[r] .= s
        features[r] .= feature_names[i]

        idx += nf[i]
    end

    DataFrame(score=scores, pval=pvals, logfc=logfc, group=groups, feature=features)
end

"""
    filter_rank_markers(de::DataFrame; pval_thresh::Real=1e-2, ngenes::Integer=typemax(Int64))

Filters and ranks a list of markers (differentially expressed genes).

**Arguments**:

    -`de`: list of markers returned by [find_markers](@ref)
    -`pval_thresh`: only keep markers with pval < pval_thresh
    -`count`: the number of highest-ranked markers to keep
    -`rankby_abs`: rank based on absolute value of the scores

**Return values**:

A `DataFrame` containing a ranked list of filtered markers.
"""
function filter_rank_markers(de::DataFrame; pval_thresh::Real=1e-2, count::Integer=typemax(Int64), rankby_abs::Bool=false)
    combine(groupby(de, :group, sort=true)) do x
        x = x[x[!,:pval] .< pval_thresh, :]
        if rankby_abs
            x[!,:abs_score] .= abs.(x[!,:score])
            sort!(x, [:abs_score, :logfc], rev=true)
            select!(x, Not(:abs_score))
        else
            sort!(x, [:score, :logfc], rev=true)
        end
        x[1:min(count, size(x,1)), :]
    end
end

function find_all_markers(X::Union{NamedCountMatrix, NamedDataMatrix}, idents::NamedVector{<:Integer}, nlabels::Int64=count_labels(idents);
        logfc_threshold::Real=0.25, min_pct::Real=0.1, min_diff_pct::Real=-Inf, only_pos::Bool=false, log::Bool=false, method=:wilcoxon)
    sel = prefilter_markers(X, idents, nlabels, logfc_threshold=logfc_threshold, min_pct=min_pct, min_diff_pct=min_diff_pct, only_pos=only_pos, log=log)
    de = find_markers(X, idents, method=method, log=log, selection=sel)
    filter_rank_markers(de)
end
