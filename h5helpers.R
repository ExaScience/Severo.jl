h5writeSparse <- function(obj, fname, dname, counts=F) {
	h5createGroup(fname, dname)
	suppressWarnings({
		x <- obj@x
		if(counts) x <- as.integer(x)
		h5write(x, fname, paste(dname, "data", sep="/"))
		h5write(obj@i, fname, paste(dname, "indices", sep="/"))
		h5write(obj@p, fname, paste(dname, "indptr", sep="/"))
		h5write(obj@Dim, fname, paste(dname, "shape", sep="/"))
		h5write(colnames(obj), fname, paste(dname, "colnames", sep="/"))
		h5write(rownames(obj), fname, paste(dname, "rownames", sep="/"))
	})
}

h5writeDataFrame <- function(obj, fname, dname) {
	h5createGroup(fname, dname)
	for(n in colnames(obj)) {
		x <- obj[,n]
		if(is.factor(x)) {
			x <- levels(x)[x]
		}

		h5write(x, fname, paste(dname, n, sep="/"))
	}
}
