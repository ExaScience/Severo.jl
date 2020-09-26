suppressPackageStartupMessages({
	library(rhdf5)
	library(SingleCellExperiment)
})

source("h5helper.R")

iname="input.rds"
oname="sce.h5"

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

sce <- readRDS(iname)
dat <- counts(sce)
x <- as(dat, "dgCMatrix")

cat("percentage nonzeros", length(x@x) / prod(x@Dim), "\n")

invisible(h5createFile(oname))
h5writeSparse(x, oname, "counts")
x <- colData(sce)
h5writeDataFrame(x, oname, "metadata")
