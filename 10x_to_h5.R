library(Matrix)
library(rhdf5)

source("h5helpers.R")

sample.dir="."
oname="10x.h5"

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

is_v3 <- file.exists(paste0(sample.dir, "/features.tsv.gz"))
suffix <- ifelse(is_v3, ".gz", "")
features.file <- ifelse(is_v3, paste0(sample.dir, "/features.tsv.gz"),
								paste0(sample.dir, "/genes.tsv")
)
matrix.file <- paste0(sample.dir, "/matrix.mtx", suffix)
barcodes.file <- paste0(sample.dir, "/barcodes.tsv", suffix)

rawdata <- readMM(matrix.file)
rawdata <- as(rawdata, "dgCMatrix")

barcodes <- readLines(barcodes.file)
colnames(rawdata) <- barcodes
features <- read.delim(features.file, header = F, stringsAsFactors = F)
rownames(rawdata) <- make.unique(features[, 2])


cat("percentage nonzeros", length(rawdata@x) / prod(rawdata@Dim), "\n")

invisible(h5createFile(oname))
h5writeSparse(rawdata, oname, "counts")
