import CSV
import GZip
import HDF5: h5open, attrs, exists

struct MMParseError <: Exception
	msg::AbstractString
end

function parseMM_header(mm::IO)
	header = readline(mm)
	tokens = split(header)
	if length(tokens) != 5
		throw(MMParseError(string("Malformed MatrixMarket header: ", header)))
	end

	if tokens[1] != "%%MatrixMarket"
		throw(ParseError(string("Not a valid MatrixMarket header:", firstline)))
	end

	(obj, format, field, symm) = map(lowercase, tokens[2:5])
	if obj != "matrix"
		throw(ParseError("Unknown MatrixMarket data type: $obj (only \"matrix\" is supported)"))
	end

	eltype = field == "real" ? Float64 :
		field == "integer" ? Int64 :
		throw(ParseError("Unsupported field $field (only real and integer are supported)"))

	if symm != "general"
		throw(ParseError("Unknown matrix symmetry: $symm (only \"general\" is supported)"))
	end

	if format != "coordinate"
		throw(ParseError("Unsupported matrix format: $format (only \"coordinate\" is supported)"))
	end

	eltype
end

function skipMM_comments(mm::IO)
	ll = readline(mm)
	while length(ll) == 0 || ll[1] == '%'
		ll = readline(mm)
	end
	ll
end

function parseMM_comments(mm::IO)
	comments = []

	ll = readline(mm)
	while length(ll) == 0 || ll[1] == '%'
		push!(comments, ll)
		ll = readline(mm)
	end

	ll, comments
end

function readMM(mm::IO; read_comments::Bool=false)
	eltype = parseMM_header(mm)

	ll, comments = if read_comments
		parseMM_comments(mm)
	else
		skipMM_comments(mm), nothing
	end

	parseint(x) = parse(Int64, x)

	# Read matrix dimensions (and number of entries)
	dd = map(parseint, split(ll))
	if length(dd) != 3
		throw(ParseError(string("Could not read in matrix dimensions: ", ll)))
	end

	rows, cols, entries = dd
	r = Vector{Int}(undef, entries)
	c = Vector{Int}(undef, entries)
	v = Vector{eltype}(undef, entries)

	for i in 1:entries
		ll = readline(mm)

		x = split(ll)
		if length(x) != 3
			throw(ParseError(string("Could not read matrix entry: ", ll)))
		end

		r[i] = parseint(x[1])
		c[i] = parseint(x[2])
		v[i] = parse(eltype, x[3])
	end

	X = sparse(r, c, v, rows, cols)
	read_comments ? (X, comments) : X
end

function readMM(fname::AbstractString; kw...)
	io = if endswith(fname, ".gz")
		GZip.open(fname, "r")
	else
		open(fname, "r")
	end

	try
		readMM(io; kw...)
	finally
		close(io)
	end
end

struct ParseError_10X <: Exception
	msg::AbstractString
end

function _read_10X_h5(fname::String, dataset::String="/mm10")
	h5open(fname, "r") do f
		feature_slot = if !exists(attrs(f), "PYTABLES_FORMAT_VERSION")
			"/features/name"
		else
			"/gene_names"
		end

		if ! exists(f, dataset)
			throw(ParseError_10X("Dataset $dataset does not exist in $fname"))
		end

		try
			p = read(f, string(dataset, "/indptr"))
			i = read(f, string(dataset, "/indices"))
			x = read(f, string(dataset, "/data"))
			dim = read(f, string(dataset, "/shape"))
			features = read(f, string(dataset, feature_slot))
			barcodes = read(f, string(dataset, "/barcodes"))
			X = SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
			X, features, barcodes
		catch e
			if isa(e, ErrorException) # probably HDF5 error
				throw(ParseError_10X("Failed to load dataset $dataset: $(e.msg)"))
			else
				rethrow(e)
			end
		end
	end
end

function readDelim(fname::AbstractString; kw...)
	io = if endswith(fname, ".gz")
		GZip.open(fname, "r")
	else
		open(fname, "r")
	end

	try
		CSV.read(io; kw...)
	finally
		close(io)
	end
end

function _read_10X(dirname::AbstractString)
	if ! isdir(dirname)
		throw(ParseError_10X("Directory $dirname does not exist"))
	end

	is_v3 = isfile(joinpath(dirname, "features.tsv.gz"))
	feature_file = joinpath(dirname, is_v3 ? "features.tsv.gz" : "genes.tsv")
	barcodes_file = joinpath(dirname, is_v3 ? "barcodes.tsv.gz" : "barcodes.tsv")
	matrix_file = joinpath(dirname, is_v3 ? "matrix.mtx.gz" : "matrix.mtx")

	if !(isfile(feature_file) && isfile(barcodes_file) && isfile(matrix_file))
		throw(ParseError_10X("Directory $dirname does not contain all components: $feature_file, $barcodes_file, $matrix_file"))
	end

	X = readMM(matrix_file)
	barcodes = readlines(barcodes_file)
	features = readDelim(feature_file, header=false)[:,2]
	X, features, barcodes
end

function read_10X(dirname::AbstractString; unique_features=true)
	X, features, barcodes = _read_10X(dirname)
	convert_data(X, features, barcodes, unique_features=unique_features)
end

function read_10X_h5(fname::String, dataset::String="/mm10"; unique_features=true)
	X, features, barcodes = _read_10X_h5(fname, dataset)
	convert_data(X, features, barcodes, unique_features=unique_features)
end

function convert_data(X::AbstractMatrix, features::AbstractVector, barcodes::AbstractVector; unique_features=true)
	if unique_features
		make_unique!(features, features)
	end

	NamedArray(copy(X'), (barcodes, features), (:cells, :features))
end
