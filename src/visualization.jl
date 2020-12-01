import Plots: plot, scatter
import StatsPlots: boxplot, violin

"""
    plot_highly_expressed_genes(X::NamedCountMatrix, n::Int64; dropfeatures::Union{Nothing, AbstractArray}=nothing)

Plot the features with the highest average expression across all cells, along with their expression in each individual cell.

**Arguments**:

    -`X`: the count matrix
    -`n`: the number of the most expressed features to show
    -`dropfeatures`: array with names, indices or bits indicating features to drop when plotting

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
    labels = names(top_counts, 2)
    boxplot(top_counts, xticks=(1:n, labels), legend=nothing)
end

"""
    plot_loadings(em::LinearEmbedding; dims::AbstractVector{<:Integer}=1:6, nfeatures::Integer=10)

Visualize top genes associated with reduction components

**Arguments**:

    -`em`: a linear embedding
    -`dims`: which components to display
    -`nfeatures`: number of genes to display

**Return value**:

A plot object
"""
function plot_loadings(em::LinearEmbedding; dims::AbstractVector{<:Integer}=1:6, nfeatures::Integer=10)
    loadings = em.basis

    featn, dimn = names(loadings)
    p = map(dims) do dim
        x = loadings.array[:,dim]
        features = partialsortperm(x, 1:nfeatures, rev=true)
        n = featn[features]
        scatter(x[features], 1:nfeatures, yticks=(1:nfeatures, n), yrotation=90, legend=nothing, title=dimn[dim])
    end
    plot(p...)
end

"""
    plot_embedding(em::LinearEmbedding)

Plots the output of a dimensional reduction technique on a 2D scatter plot where each point is a
cell and it's positioned based on the cell embeddings determined by the reduction technique.

**Arguments**:

    -`em`: a linear embedding

**Return value**:

A plot object

"""
function plot_embedding(em::LinearEmbedding)
    labels = names(em.coordinates,2)
    scatter(em.coordinates[:,1], em.coordinates[:,2], xlabel=labels[1], ylabel=labels[2], legend=nothing)
end

"""
    plot_elbow(em::LinearEmbedding)

Plots the standard deviations of the principle components for easy identification of an elbow in the graph.

**Arguments**:

    -`em`: a linear embedding

**Return value**:

A plot object

"""
function plot_elbow(em::LinearEmbedding; screeplot::Bool=true)
    if screeplot
        x = cumsum(em.stdev) ./ sum(em.stdev)
        scatter(x, xlabel="Component", ylabel="Explained variability", legend=nothing)
    else
        scatter(em.stdev, xlabel="Component", ylabel="Standard deviation", legend=nothing)
    end
end

function plot_rank_features_group(dx::DataFrame, group::Integer; features::Union{Nothing,AbstractArray}, ngenes::Integer=20)
    if features === nothing
        features = 1:ngenes
    end


end

function plot_rank_features_group(dx::DataFrame, groups::AbstractVector{<:Integer}; features::Union{Nothing,AbstractArray}, ngenes::Integer=20)
end

function plot_violin(X::NamedCountMatrix, feature::Union{Integer, AbstractString}, lbls::AbstractVector{<:Integer})
    x = X[:,feature]
    feature_name = if feature isa AbstractString
        feature
    else
        names(X,2)[feature]
    end

    violin(lbls, x, legend=nothing, ylabel=feature_name)
end

function plot_violin(X::NamedCountMatrix, features::AbstractArray, lbls::AbstractVector{<:Integer})
    p = map(features) do feat
        plot_violin(X, feat, lbls)
    end

    plot(p...)
end

export plot_highest_expressed, plot_embedding, plot_loadings, plot_elbow, plot_violin
