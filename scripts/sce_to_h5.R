suppressPackageStartupMessages({
	library(rhdf5)
	library(SingleCellExperiment)
})

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

iname="input.rds"
oname <- NULL

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

if(is.null(oname)) {
	oname <- paste0(tools::file_path_sans_ext(iname), ".h5")
}

sce <- readRDS(iname)
dat <- counts(sce)
x <- as(dat, "dgCMatrix")

cat("percentage nonzeros", length(x@x) / prod(x@Dim), "\n")

invisible(h5createFile(oname))
h5writeSparse(x, oname, "counts", counts=T)
x <- colData(sce)
h5writeDataFrame(x, oname, "metadata")
