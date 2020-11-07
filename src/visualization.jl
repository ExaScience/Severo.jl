"""
    plot_highly_expressed_genes(X::NamedCountMatrix, n::Int64; dropfeatures::Union{Nothing, AbstractArray}=nothing)

Plot the features with the highest average expression across all cells, along with their expression in each individual cell.

**Arguments**:

    -``X``: the count matrix
    -``n``: the number of the most expressed features to show
    -``dropfeatures``: array with names, indices or bits indicating features to drop when plotting

**Return value**:

A plot object
"""
function plot_highest_expressed(X::NamedCountMatrix; n::Int64=10, dropfeatures::Union{Nothing, AbstractArray}=nothing)
    norm = normalize(X, scale_factor=100, method=:relativecounts)
    if dropfeatures !== nothing
        norm = norm[:, Not(dropfeatures)]
    end

    mean_percent = vec(mean(norm, dims=1))
    top_idx = partialsortperm(mean_percent, 1:n, rev=true)

    top_counts = norm[:, top_idx]
    labels = reshape(names(top_counts, 2), 1, :)
    boxplot(top_counts, labels=labels)
end
