import Loess: loess, predict
import Statistics: var

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
