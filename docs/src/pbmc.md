# Preprocessing and clustering 3k PBMCs

This showcase reproduces Seurat's [Guided Clustering Tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html)

## Overview

```@contents
Pages = ["pbmc.md"]
```

## Loading data

The data consists of 3k PBMCs from a Healthy Donor and is freely available from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k).

You can either download the data manually
```bash
mkdir pbmc3k
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
cd pbmc3k; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz --strip-components 2
```

The `read_10X` function reads the data from the cellranger 10X pipeline and returns a labelled count matrix. Each entry indicated the number of molecules
detected for each feature/gene (columns) and each cell (rows).

```julia
using Cell
X = read_10X("pbmc3k/")
```

!!! warning
    The meaning of the rows and columns is different from the representation used by other packages like [Seurat](https://satijalab.org/seurat/)

Alternatively, the ```dataset``` function can be used to load a dataset from a predefined collection. For example, the PBMC collection is known by `Cell.jl`
and can be easily loaded as follows:

```@example pbmc
using Cell
X = dataset("PBMC", "3k")

# The matrix can be indexed using names or indices. For instance,
# we can look at specific genes in the first thirty cells
X[1:30, ["CD3D", "TCL1A", "MS4A1"]]
```

The count data is stored in a sparse matrix format, only storing non-zero elements of the matrix. Any values not shown are zero.

```@docs
read_data
read_10X
dataset
```

## Preprocessing

```@example pbmc
using Plots
using StatsPlots
pyplot() # hide
plot_highest_expressed(X)
savefig("he-plot.svg"); nothing # hide
```

![](he-plot.svg)


### Filtering

Two filtering function are available for filtering cells and features/genes based on basic criteria such as the number of features
detected, the number of cells for which a feature is detected and the total number of counts.

```@example pbmc
X = filter_cells(X, min_features=200)
X = filter_features(X, min_cells=3)
```

This filters out cells with less than 200 features/genes and features detected in less than 3 cells. For convenience, `filter_counts`
combines these two into a single function, similar to Seurat. Beware that the order in which you call these functions can matter.

```@docs
filter_features
filter_cells
filter_counts
```

### Normalization

After removing unwanted cells from the dataset, the next step is to normalize the data. The basic method normalizes the feature expression
measurements for each cell by the total expression, multiplies this by a scale factor and optionally log-transforming the result.

```@example pbmc
Y = normalize(X, method=:lognormalize, scale_factor=1e4)
```

```@docs
normalize
```

### Identifying highly variable features

Next, we select a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).

```@example pbmc
hvf = find_variable_features(X, 2000; method=:vst)
```

```@docs
find_variable_features
```

## Dimensionality reduction

### Scaling

Prior to dimensional reduction techniques like PCA, it's a good idea to `scale` the data. Scaling performs two basic operations:

- Shifts the expression of each gene, so that the mean expression across cells is 0
- Scales the expression of each gene, so that the standard deviation across cells is 1. This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
- Clips values exceeding standard deviation of `scale_max`

```@example pbmc
Y = Y[:,hvf] # only use highly-variable features
S = scale(Y; scale_max=10)
```

```@docs
scale
```

### Principal component analysis

Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and de-noises the data.

```@example pbmc
em = embedding(S, 15, method=:pca)
```

```@example pbmc
plot_elbow(em)
savefig("elbow-plot.svg"); nothing # hide
```

![](elbow-plot.svg)

```@docs
embedding
```

## Clustering the cells

Briefly, the clustering embed cells in a graph structure with edges between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'.

### Computing the neighborhood graph

First, we need to find the k-nearest neighbors of all the cells based on the euclidean distance in PCA space, and then compute a graph based on the shared overlap in the local neighborhoods between any two cells (Jaccard similarity). These step are performed using the `nearest_neighbours` and `jaccard_index` functions. For convenience, the two steps are also combined in the `shared_nearest_neighours` function.

```@example pbmc
nn = nearest_neighbours(em, 20, 1:10)
snn = jaccard_index(nn, prune=1/15)
# or as a single step
#snn = nearest_neighbours(em, 20, 1:10, prune=1/15)
```

```@docs
nearest_neighbours
jaccard_index
shared_nearest_neighbours
```

### Clustering the neighborhood graph

Clustering is performed using a graph-based method called `Louvain clustering` which tries to detect communities based on optimizing a modularity function.
Since this is a greedy algorithm, which can get stuck in local minima based on random initialization, multiple restarts are required to find a more "global" minima.

```@example pbmc
lbls = cluster(snn)
```

```@docs
cluster
```

## Finding differentially expressed features

We can find markers that define clusters via differential expression. It identifies positive and negative markers for a single cluster compared to the other cells.
`find_markers` performs this process for every cluster. The function returns a `DataFrame` containing for each marker: the p-value, the log-foldchange and statistical score.

```@example pbmc
dx = find_markers(X, lbls; method=:wilcoxon)
dx = filter_rank_markers(dx)
dx[1:10, :]
```

```@docs
find_markers
```

```@example pbmc
plot_violin(X, ["CST3", "NKG7", "PPBP", "S100A4"], lbls)
savefig("violin-plot.svg"); nothing # hide
```

![](violin-plot.svg)

