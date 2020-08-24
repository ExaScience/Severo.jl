library(proxyC)
library(irlba)
library(Matrix)
library(rhdf5)
source("tictoc.R")

fname <- "/data/thaber/1M_neurons_filtered_gene_bc_matrices_cut.h5"

load_h5 <- function(fname) {
	A <- new("dgCMatrix")
	A@i <- as.integer(h5read(fname, "mm10/indices"))
	A@p <- as.integer(h5read(fname, "mm10/indptr"))
	A@x <- as.numeric(h5read(fname, "mm10/data"))
	A@Dim <- as.integer(h5read(fname, "mm10/shape"))
	rownames(A) <- h5read(fname, "mm10/gene_names")
	colnames(A) <- h5read(fname, "mm10/barcodes")
	A <- t(A)
}

A <- load_h5(fname)
elms <- h5ls(fname, recursive=F)$name

tic()
mu <- colMeans(A)
std <- colSds(A)
cat("mean/var", toc(F, T), "\n")

if( "features" %in% elms ) {
	features <- h5read(fname, "features")
	i <- match(features, colnames(A))
} else {
	i <- sort(sample(which(std > 0), 2000))
}

B <- A[,i]

tic()
pca.res <- irlba(B, 100, center=mu[i], scale=std[i])
cat("pca irlba sparse", toc(F, T), "\n")

if( "scaled" %in% elms ) {
	S <- h5read(fname, "scaled")
	S <- t(S)
	rownames(S) <- rownames(B)
	colnames(S) <- colnames(B)
} else {
	tic()
	S <- sweep(sweep(B, 2, mu[i], FUN=`-`), 2, std[i], FUN=`/`)
	cat("create dense matrix", toc(F, T), "\n")
}

tic()
pca.res2 <- irlba(S, 100)
cat("pca irlba dense", toc(F, T), "\n")

tic()
prc.res <- prcomp(S, rank.=100); toc()
cat("pca prcomp dense", toc(F, T), "\n")

tic()
G = crossprod(S, S)
pca.gram <- irlba(G, 100)
U <- S %*% pca.gram$v %*% diag(1 / pca.gram$d)
cat("pca gram", toc(F, T), "\n")

