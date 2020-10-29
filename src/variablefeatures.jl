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

function variance_stabilizing_transformation(A::SparseMatrixCSC{<:Integer}; loess_span=0.5)
    mu, std = mean_std(A)
    non_const = std .> 0

    xs, ys = log10.(mu[non_const]), log10.(std[non_const])
    model = loess(xs, ys, span=loess_span)

    expected_std = std # reuse std memory
    expected_std[non_const] = 10 .^ predict(model, xs)

  # workaround for bug in loess
    expected_std[isnan.(expected_std)] .= 0.0

    standardized_var_clipped(A, mu, expected_std)
end

qnorm(q::Real) = quantile(Normal(), q)

# In the first stage, digital gene expression matrices were column-normalized. Cells with fewer than
# 400 expressed genes were removed from analysis. To identify a set of highly variable genes,
# we first calculated the average mean and variance of each gene, and selected genes that
# were: (1) 0.1 log10 units above the expected variance for a perfectly Poisson-distributed
# gene of equivalent mean expression; and (2) above a Bonferroni-corrected 99% confidence
# interval defined by a normal approximation of a Poisson distribution. Saunders et al. 2018
function select_features_saunders(counts::SparseMatrixCSC{<:Integer}, norm::SparseMatrixCSC{<:Real}; var_thresh=0.1, alpha_thresh=0.99)
    ncells, ngenes = size(counts)
    trx_per_cell = vec(sum(counts, dims=2))
    gene_expr_mean, gene_expr_var = mean_var(norm)
    nolan_constant = mean(1 ./ trx_per_cell)
    alphathresh_corrected = alpha_thresh / ngenes
    genemeanupper = gene_expr_mean .+ qnorm(1 - alphathresh_corrected / 2) * sqrt.(gene_expr_mean * nolan_constant / ncells)
    basegenelower = log10.(gene_expr_mean .* nolan_constant)

    #findall(((gene_expr_var ./ nolan_constant) .> genemeanupper) .& (log10.(gene_expr_var) .> basegenelower .+ var_thresh))
    metric = zeros(ngenes)
    J = ((gene_expr_var ./ nolan_constant) .> genemeanupper)
    metric[J] = log10.(gene_expr_var[J]) .- basegenelower[J]
    metric
end

function find_variable_features(counts::NamedCountMatrix; method=:vst, nfeatures=2000, kw...)
    if isa(method, AbstractString)
        method = Symbol(method)
    end

    metric = if method == :vst
        variance_stabilizing_transformation(counts.array, kw...)
    elseif method == :saunders
        norm = if :norm in keys(kw)
            norm = kw[:norm]
            isa(norm, NamedArray) ? norm.array : norm
        else
            row_norm(counts.array)
        end

        var_thresh = get(kw, :var_thresh, 0.1)
        alpha_thresh = get(kw, :alpha_thresh, 0.1)
        select_features_saunders(counts.array, norm; alpha_thresh=alpha_thresh, var_thresh=var_thresh)
    else
        error("unknown selection method: $method")
    end

    selected = partialsortperm(metric, 1:nfeatures, rev=true)
    NamedArray(selected, (names(counts,2)[selected],), (dimnames(counts,2),))
end
