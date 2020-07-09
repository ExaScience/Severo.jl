library(dplyr)
library(future)
library(Seurat)
source('tictoc.R')
#pdf('bio.pdf')

#options(future.globals.maxSize=2*1024*1024*1024)
#plan("multiprocess", workers=4)

tic()
pbmc.data <- Read10X(data.dir=".")
cat("load data", toc(F, T), "\n")

pbmc <- CreateSeuratObject(counts = pbmc.data, project="pbmc3k", min.cells=3, min.features=200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
cat("CreateSeuratObject", toc(F, T), "\n")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
cat("NormalizeData", toc(F, T), "\n")

pbmc <- FindVariableFeatures(pbmc, selection.method="vst", nfeatures=2000)
cat("FindVariableFeatures", toc(F, T), "\n")

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features=all.genes)
cat("ScaleData", toc(F, T), "\n")

pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, rev.pca=T)
cat("RunPCA", toc(F, T), "\n")

#ElbowPlot(pbmc, ndims = 100)
#DimHeatmap(pbmc, dims = c(1:3, 70:75), cells = 500, balanced = TRUE)

tic()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
cat("FindNeighbors", toc(F, T), "\n")
pbmc <- FindClusters(pbmc, resolution = 0.5)
cat("FindClusters", toc(F, T), "\n")

pbmc <- FindAllMarkers(pbmc, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
cat("FindAllMarkers", toc(F, T), "\n")

#pbmc <- RunUMAP(pbmc, dims = 1:75, min.dist = 0.75)
#cat("RunUMAP", toc(F, T), "\n")

#library(ggplot2)
#p1 <- DimPlot(pbmc, reduction = "tsne", pt.size = 0.1) + ggtitle(label = "FIt-SNE")
##p2 <- DimPlot(pbmc, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP")
#p1 <- AugmentPlot(plot = p1)
#p2 <- AugmentPlot(plot = p2)
##CombinePlots(plots = list(p1, p2), legend = "none")
#CombinePlots(plots = p1)

#p3 <- FeaturePlot(pbmc, features = c("S100a9", "Sftpc"), reduction = "tsne", pt.size = 0.1, combine = FALSE)
#p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
#CombinePlots(plots = p3)

#dev.off()
