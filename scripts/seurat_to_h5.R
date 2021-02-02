suppressPackageStartupMessages({
	library(rhdf5)
	library(Seurat)
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

obj <- readRDS(iname)
x <- GetAssayData(obj, "counts")
meta.data <- obj[[]]
id <- Idents(obj)
hvf <- VariableFeatures(obj)

invisible(h5createFile(oname))
h5writeSparse(x, oname, "counts", counts=T)
h5writeDataFrame(meta.data, oname, "metadata")

if(!is.null(id)) {
	if(is.factor(id)) {
		id <- levels(id)[id]
	}
	h5write(id, oname, "idents")
}

if(!is.null(hvf)) {
	h5write(hvf, oname, "hvf")
}
