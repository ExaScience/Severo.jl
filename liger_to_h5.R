library(liger)
library(rhdf5)

source("h5helpers.R")

iname="input.rds"
oname=NULL
store_intermediates=TRUE
store_nmf=TRUE
store_clusters=TRUE

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
	eval(parse(text = args[[i]]))
}

if(is.null(oname)) {
	oname <- paste0(tools::file_path_sans_ext(iname), ".h5")
}

liger <- readRDS(iname)
raw.data <- liger@raw.data
datasets <- names(raw.data)

invisible(h5createFile(oname))

for(ds in datasets) {
	h5createGroup(oname, ds)

	counts <- raw.data[[ds]]
	h5writeSparse(counts, oname, paste0(ds, "/counts"), counts=T)

	if(store_intermediates) {
		h5writeSparse(liger@norm.data[[ds]], oname, paste0(ds, "/norm"))
		h5write(liger@scale.data[[ds]], oname, paste0(ds, "/scale"))
	}
}

if(store_intermediates) {
	h5write(liger@var.genes, oname, "/features")
}

if(store_nmf) {
	h5createGroup(oname, "nmf")
	h5write(liger@W, oname, "nmf/W")

	h5createGroup(oname, "nmf/H")
	h5createGroup(oname, "nmf/V")
	for(ds in datasets) {
		h5write(liger@H[[ds]], oname, paste0("nmf/H/", ds))
		h5write(liger@V[[ds]], oname, paste0("nmf/V/", ds))
	}
}

if(store_clusters) {
	h5write(as.integer(liger@clusters), oname, "idents")
}
