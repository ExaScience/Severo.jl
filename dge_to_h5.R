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

iname <- "input.raw.dge.txt.gz"
oname <- NULL

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

strEndsWith <- function(theString,theExt) {
    return(substring(theString,1+nchar(theString)-nchar(theExt))==theExt)
}

strStartsWith <- function(theString, thePrefix) {
    return(strtrim(theString, nchar(thePrefix)) == thePrefix)
}

removeSuffix <- function(str, suffix) {
    if (!strEndsWith(str, suffix)) {
        stop(paste0("'", str, "' does not end with '", suffix, "'"))
    }
    return(strtrim(str, nchar(str) - nchar(suffix)))
}

maybeRemoveSuffix <- function(str, suffix) {
    if (!strEndsWith(str, suffix)) {
        return(str)
    } else {
        return(strtrim(str, nchar(str) - nchar(suffix)))
    }
}

#' returns list(genes, cell_barcodes)
loadSparseDgeNames <- function(file) {
	conn <- file(file, "r")
	genes <- c()
	cell_barcodes <- c()
	while(TRUE) {
		line <- readLines(con=conn, n=1, ok=FALSE)
		if (!strStartsWith(line, '%')) {
			break
		}
		if (strStartsWith(line, '%%GENES\t')) {
			these_genes <- strsplit(line, "\t", fixed=TRUE)[[1]]
			genes <- append(genes, these_genes[2:length(these_genes)])
		} else if (strStartsWith(line, '%%CELL_BARCODES\t')) {
			these_cells <- strsplit(line, "\t", fixed=TRUE)[[1]]
			cell_barcodes <- append(cell_barcodes, these_cells[2:length(these_cells)])
		}
	}
	close(conn)
	return(list(genes=genes, cell_barcodes=cell_barcodes))
}

#' Load a sparse DGE (either real or integer), that may have gene names and cell_barcode names in header
#' @return A sparse matrix in triplet form (dgTMatrix)
loadSparseDge <- function(file) {
	genes_and_cell_barcodes <- loadSparseDgeNames(file)
	ret <- readMM(file)
	rownames(ret) <- genes_and_cell_barcodes$genes
	colnames(ret) <- genes_and_cell_barcodes$cell_barcodes
	return(ret)
}

bname <- removeSuffix(maybeRemoveSuffix(iname, ".gz"), ".raw.dge.txt")
if(is.null(oname)) {
	oname <- paste0(tools::file_path_sans_ext(bname), ".h5")
}

x <- loadSparseDge(iname)
x <- as(x, "dgCMatrix")

cat("percentage nonzeros", length(x@x) / prod(x@Dim), "\n")

invisible(h5createFile(oname))
h5writeSparse(x, oname, "counts")

cluster.file <- paste0(bname, ".cluster.assign.RDS")
if( file.exists(cluster.file) ) {
	cl <- readRDS(cluster.file)
	h5write(as.integer(cl), oname, "/idents")
}
