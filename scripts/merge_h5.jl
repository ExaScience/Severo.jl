# copyright imec - evaluation license - not for distribution

import SparseArrays: SparseMatrixCSC, getcolptr, nonzeros, rowvals
import HDF5: h5read, h5write
import Printf: @sprintf

read_gene_names(fname, dataset="/counts") = h5read(fname, @sprintf("%s/rownames", dataset))
read_barcodes(fname, dataset="/counts") = h5read(fname, @sprintf("%s/colnames", dataset))

function read_sparse(fname, dataset="/counts")
	p = h5read(fname, @sprintf("%s/indptr", dataset))
	i   = h5read(fname, @sprintf("%s/indices", dataset))
	x = h5read(fname, @sprintf( "%s/data", dataset))
	dim = h5read(fname, @sprintf("%s/shape", dataset))
	SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
end

function read_sparse_named(fname, dataset="/counts")
	A = read_sparse(fname, dataset)
	rownames = h5read(fname, @sprintf("%s/rownames", dataset))
	colnames = h5read(fname, @sprintf("%s/colnames", dataset))
	A, rownames, colnames
end

function write_sparse(fname, X::SparseMatrixCSC, dataset="/counts")
	h5write(fname, @sprintf("%s/indptr", dataset), getcolptr(X) .- 1)
	h5write(fname, @sprintf("%s/indices", dataset), rowvals(X) .- 1)
	h5write(fname, @sprintf( "%s/data", dataset), nonzeros(X))
	h5write(fname, @sprintf("%s/shape", dataset), collect(size(X)))
end

function write_sparse_named(fname, X, rownames, colnames, dataset="/counts")
	write_sparse(fname, X, dataset)
	h5write(fname, @sprintf("%s/rownames", dataset), rownames)
	h5write(fname, @sprintf("%s/colnames", dataset), colnames)
end

function union_index(X::AbstractArray)
	x = foldl(union, X)
	d = Dict{eltype(x), eltype(inds)}(zip(x, keys(x)))
	x, d
end

function merge_datasets(files, dataset="/counts"; min_genes=nothing)
	# calculate union of all genes
	all_genes = map(read_gene_names, files)
	common_genes, index = union_index(all_genes)
	Xs,barcodes = zip(map(files) do fname
		X, genes, barcodes = read_sparse_named(fname, dataset)
		inds = [get(index, g, 0) for g in genes]

		if min_genes !== nothing
			CI = vec(sum(X, dims=1)) .> min_genes
			X[inds, CI], barcodes[CI]
		else
			X[inds, :], barcodes
		end
	end...)
	X = hcat(Xs...)
	barcodes = vcat(barcodes...)
	X, common_genes, barcodes
end

function memory_usage()
	mem = open("/proc/self/statm", "r") do io
		parse.(Int64, split(readline(io))) .* 4096
	end
	mem[2] / 1024/1024/1024
end

oname = ARGS[1]
files = ARGS[2:end]
# files = joinpath.("/data/thaber/TabulaMuris/ss2/", filter(x -> endswith(x, ".h5"), readdir("/data/thaber/TabulaMuris/ss2/")))

X, genes, barcodes = merge_datasets(files)
write_sparse_named(oname, X, genes, barcodes)
