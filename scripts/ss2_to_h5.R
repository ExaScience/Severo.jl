library(rhdf5)
library(Matrix)

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

iname="input.csv"
oname=NULL

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

if(is.null(oname)) {
	oname <- paste0(tools::file_path_sans_ext(iname), ".h5")
}

X <- read.csv(iname, row.names=1)
counts <- as(as.matrix(X), "dgCMatrix")

invisible(h5createFile(oname))
h5writeSparse(counts, oname, "/counts", counts=T)

