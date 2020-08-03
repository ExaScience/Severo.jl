library(proxyC)
library(irlba)
library(Matrix)
library(rhdf5)
source("tictoc.R")

fname <- "/data/thaber/1M_neurons_filtered_gene_bc_matrices_h5.h5"

load_h5 <- function(fname) {
	A <- new("dgCMatrix")
	A@i <- as.integer(h5read(fname, "mm10/indices"))
	A@p <- as.integer(h5read(fname, "mm10/indptr"))
	A@x <- as.numeric(h5read(fname, "mm10/data"))
	A@Dim <- as.integer(h5read(fname, "mm10/shape"))
	rownames(A) <- h5read(fname, "mm10/genes")
	colnames(A) <- h5read(fname, "mm10/barcodes")
	A <- t(A)
}

A <- load_h5(fname)
elms <- h5ls(fname, recursive=F)$name

if( "features" %in% elms ) {
	features <- h5read(fname, "features")
	i <- match(features, colnames(A))
	B <- A[,i]
} else {
	B <- A
	i <- 1:ncol(A)
}

mu <- colMeans(A)
std <- colSds(A)

tic(); pca.res <- irlba(B, 100, center=mu[i], scale=std[i]); toc()

if( "scaled" %in% elms ) {
	S <- h5read("/data/thaber/pca.h5", "scaled")
	S <- t(S)
	rownames(S) <- rownames(B)
	colnames(S) <- colnames(B)
} else {
	S <- sweep(sweep(B, 2, mu[i], FUN=`-`), 2, std[i], FUN=`/`)
}

tic(); pca.res2 <- irlba(S, 100); toc()
tic(); prc.res <- prcomp(S, rank.=100); toc()
