import Loess: loess, predict
import Statistics: var, mean
import Distributions: quantile, Normal

standardize_clip(x, mu, std, vmax) = min((x - mu) / std, vmax)

function standardized_var_clipped(A::SparseMatrixCSC{T}, mu::Vector{R}, std::Vector{R}; vmax=sqrt(size(A,1))) where {T <: Signed, R <: Real}
	var_stand = zeros(size(A,2))
	for (i,x) in enumerate(eachcol(A))
		std[i] == 0.0 && continue
		var_stand[i] = (sum(standardize_clip.(nonzeros(x), mu[i], std[i], vmax).^2) + (length(x) - nnz(x)) * standardize_clip(0, mu[i], std[i], vmax)^2) / (length(x)-1)
	end
	var_stand
end

function find_variable_features(A::SparseMatrixCSC{T}, nfeatures=2000; loess_span=0.5) where {T <: Signed}
	mu, std = mean_std(A)
	non_const = std .> 0

	xs, ys = log10.(mu[non_const]), log10.(std[non_const])
	model = loess(xs, ys, span=loess_span)

	expected_std = std # reuse std memory
	expected_std[non_const] = 10 .^ predict(model, xs)

  # workaround for bug in loess
	expected_std[isnan.(expected_std)] .= 0.0

	var_std = standardized_var_clipped(A, mu, expected_std)
	partialsortperm(var_std, 1:nfeatures, rev=true)
end

qnorm(q::Real) = quantile(Normal(), q)

function select_features_saunders(counts::SparseMatrixCSC{T}, norm::SparseMatrixCSC{R}; var_thresh=0.1, alpha_thresh=0.99) where {T <: Signed, R <: Real}
	ncells, ngenes = size(counts)
	trx_per_cell = vec(sum(counts, dims=2))
	gene_expr_mean, gene_expr_var = mean_var(norm)
	nolan_constant = mean(1 ./ trx_per_cell)
	alphathresh_corrected = alpha_thresh / ngenes
	genemeanupper = gene_expr_mean .+ qnorm(1 - alphathresh_corrected / 2) * sqrt.(gene_expr_mean * nolan_constant / ncells)
	basegenelower = log10.(gene_expr_mean .* nolan_constant)

	findall(((gene_expr_var ./ nolan_constant) .> genemeanupper) .& (log10.(gene_expr_var) .> basegenelower .+ var_thresh))
end
