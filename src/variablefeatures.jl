# copyright imec - evaluation license - not for distribution

import Loess: loess, predict
import Statistics: var, mean
import Distributions: quantile, Normal

standardize_clip(x, mu, std, vmax) = min((x - mu) / std, vmax)

function standardized_var_clipped(A::SparseMatrixCSC{<:Integer}, mu::Vector{R}, std::Vector{R}; vmax=sqrt(size(A,1))) where {R <: Real}
    var_stand = zeros(size(A,2))
    for (i,x) in enumerate(eachcol(A))
        std[i] == 0.0 && continue
        var_stand[i] = (sum(standardize_clip.(nonzeros(x), mu[i], std[i], vmax).^2) + (length(x) - nnz(x)) * standardize_clip(0, mu[i], std[i], vmax)^2) / (length(x)-1)
    end
    var_stand
end

function nan2zero!(v::AbstractVector{T}) where {T <: AbstractFloat}
    v .= ifelse.(isnan.(v), 0.0, v)
end

function variance_stabilizing_transformation(A::SparseMatrixCSC{<:Integer}; loess_span::Real=0.5, dtype::Type{<:AbstractFloat}=Float64)
    loess_span = convert(dtype, loess_span)

    mu, std = mean_std(dtype, A)
    non_const = std .> 0

    xs, ys = log10.(mu[non_const]), log10.(std[non_const])
    model = loess(xs, ys, span=loess_span)

    expected_std = std # reuse std memory
    expected_std[non_const] = 10 .^ predict(model, xs)

  # workaround for bug in loess
    expected_std = nan2zero!(expected_std)

    standardized_var_clipped(A, mu, expected_std)
end

qnorm(q::Real) = quantile(Normal(), q)

# In the first stage, digital gene expression matrices were column-normalized. Cells with fewer than
# 400 expressed genes were removed from analysis. To identify a set of highly variable genes,
# we first calculated the average mean and variance of each gene, and selected genes that
# were: (1) 0.1 log10 units above the expected variance for a perfectly Poisson-distributed
# gene of equivalent mean expression; and (2) above a Bonferroni-corrected 99% confidence
# interval defined by a normal approximation of a Poisson distribution. Saunders et al. 2018
function select_features_saunders(counts::SparseMatrixCSC{<:Integer}, norm::SparseMatrixCSC{T}; alpha_thresh::Real=0.99) where {T<:Real}
    ncells, ngenes = size(counts)
    trx_per_cell = vec(sum(counts, dims=2))
    gene_expr_mean, gene_expr_var = mean_var(norm)
    nolan_constant = mean(1 ./ trx_per_cell)
    alphathresh_corrected = convert(T, alpha_thresh) / ngenes
    genemeanupper = gene_expr_mean .+ qnorm(1 - alphathresh_corrected / 2) * sqrt.(gene_expr_mean * nolan_constant / ncells)
    basegenelower = log10.(gene_expr_mean .* nolan_constant)

    metric = zeros(ngenes)
    J = ((gene_expr_var ./ nolan_constant) .> genemeanupper)
    metric[J] = log10.(gene_expr_var[J]) .- basegenelower[J]
    metric
end

function select_dispersion(norm::SparseMatrixCSC{<:Real}; num_bins=20, binning_method=:width)
    mu, disp = log_VMR(norm)
    nan2zero!(disp)
end

function select_meanvarplot(norm::SparseMatrixCSC{<:Real}; num_bins=20, binning_method=:width)
    mu, disp = log_VMR(norm)
    mu = nan2zero!(mu)
    disp = nan2zero!(disp)

    _, bins = cut(mu, num_bins; method=binning_method)
    bin_mean, bin_std = mean_std(disp, bins, num_bins)

    @inbounds @simd for i = 1:length(disp)
        k = bins[i]
        disp[i] = (disp[i] - bin_mean[k]) / bin_std[k]
    end

    nan2zero!(disp)
end

"""
    find_variable_features(counts::NamedCountMatrix, nfeatures=2000; method=:vst, kw...)

Identification of highly variable features: find features that exhibit high cell-to-cell variation in the dataset
(i.e, they are highly expressed in some cells, and lowly expressed in others).

**Arguments**:

    -`counts`: count matrix
    -`nfeatures`: the number of top-ranking features to return
    -`method`: how to choose top variable features
    -`kw`: additional keyword arguments to pass along to the method

*Methods**:

    -`:vst`: fits a line to the log(mean) - log(variance) relationship, then standardizes the features values
        using the observed mean and expected variance. Finally, feature variance is calculated using the standardized values.

        - `loess_span`: span parameter for loess regression when fitting the mean-variance relationship

    -`dispersion`: selects the genes with the highest dispersion values

    -`meanvarplot`: calculates the feature mean and dispersion, bins the mean according into `num_bins` bins.
        Finally, returns the z-scores for dispersion within each bin.

        - `num_bins`: Total number of bins to use
        - `binning_method`: Specifies how the bins should be computed. Available: `:width` for equal width and `:frequency` for equal frequency binning

**Return value**:

The `nfeatures` top-ranked features
"""
function find_variable_features(counts::NamedCountMatrix, nfeatures=2000; method=:vst, dtype::Type{<:AbstractFloat}=Float64, kw...)
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    metric = if method == :vst
        variance_stabilizing_transformation(counts.array, dtype=dtype, kw...)
    else
        norm = if :norm in keys(kw)
            norm = kw[:norm]
            isa(norm, NamedArray) ? norm.array : norm
        else
            row_norm(counts.array, one(dtype))
        end

        if method == :saunders
            alpha_thresh = get(kw, :alpha_thresh, 0.1)
            select_features_saunders(counts.array, norm; alpha_thresh=alpha_thresh)
        elseif method == :dispersion
            num_bins = get(kw, :num_bins, 20)
            binning_method = get(kw, :binning_method, :width)
            select_dispersion(norm; num_bins=num_bins, binning_method=binning_method)
        elseif method == :meanvarplot
            num_bins = get(kw, :num_bins, 20)
            binning_method = get(kw, :binning_method, :width)
            select_meanvarplot(norm; num_bins=num_bins, binning_method=binning_method)
        else
            error("unknown selection method: $method")
        end
    end

    selected = partialsortperm(metric, 1:nfeatures, rev=true)
    NamedArray(selected, (names(counts,2)[selected],), (dimnames(counts,2),))
end
