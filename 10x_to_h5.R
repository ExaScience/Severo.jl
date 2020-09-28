library(Matrix)
library(rhdf5)

thisFile <- function() {
	cmdArgs <- commandArgs(trailingOnly = FALSE)
	needle <- "--file="
	match <- grep(needle, cmdArgs)
	if (length(match) > 0) {
		# Rscript
		return(normalizePath(sub(needle, "", cmdArgs[match])))
	} else {
		# 'source'd via R console
		return(normalizePath(sys.frames()[[1]]$ofile))
	}
}

source(file.path(dirname(thisFile()), "h5helpers.R"))

folder <- NULL
sample.dir <- NULL
oname <- NULL

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

if(is.null(oname) && !is.null(sample.dir)) {
	oname <- file.path(sample.dir, "10x.h5")
}

is_10x <- function(sample.dir) {
	is_v3 <- file.exists(file.path(sample.dir, "features.tsv.gz"))
	suffix <- ifelse(is_v3, ".gz", "")
	matrix.file <- file.path(sample.dir, paste0("matrix.mtx", suffix))
	file.exists(matrix.file)
}

convert_10x <- function(sample.dir, oname, library.name=NULL) {
	cat("convert", sample.dir, oname, library.name, "\n")
	is_v3 <- file.exists(file.path(sample.dir, "features.tsv.gz"))
	suffix <- ifelse(is_v3, ".gz", "")
	features.file <- ifelse(is_v3, file.path(sample.dir, "features.tsv.gz"),
									file.path(sample.dir, "genes.tsv")
	)
	matrix.file <- file.path(sample.dir, paste0("matrix.mtx", suffix))
	barcodes.file <- file.path(sample.dir, paste0("barcodes.tsv", suffix))

	rawdata <- readMM(matrix.file)
	rawdata <- as(rawdata, "dgCMatrix")

	barcodes <- readLines(barcodes.file)
	if(!is.null(library.name)) {
		barcodes <- paste0(library.name, "_", barcodes)
	}
	colnames(rawdata) <- barcodes
	features <- read.delim(features.file, header = F, stringsAsFactors = F)
	rownames(rawdata) <- make.unique(features[, 2])


	cat("percentage nonzeros", length(rawdata@x) / prod(rawdata@Dim), "\n")

	invisible(h5createFile(oname))
	h5writeSparse(rawdata, oname, "counts", counts=T)
}

library_name_10x <- function(name) strsplit(name, "-")[[1]][2]

if(is.null(folder)) {
	convert_10x(sample.dir, oname)
} else {
	sample.dirs <- list.dirs(folder)
	for(sample.dir in sample.dirs) {
		if(is_10x(sample.dir)) {
			name <- basename(sample.dir)
			library.name <- library_name_10x(name)
			oname <- file.path(folder, paste0(tolower(name), ".h5"))
			convert_10x(sample.dir, oname, library.name=library.name)
		}
	}
}
