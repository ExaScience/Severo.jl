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

"""
    read_10X(dirname::AbstractString; unique_features=true)

Read count matrix from 10X genomics

**Arguments**:

- `dirname`: path to directory containing matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv from 10X
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_10X(dirname::AbstractString; unique_features::Bool=true)
    X, features, barcodes = _read_10X(dirname)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    read_10X_h5(fname::String, dataset::String="/mm10"; unique_features=true)

Read count matrix from 10X CellRanger hdf5 file.

**Arguments**:

- `fname`: path to hdf5 file
- `dataset`: name of dataset to load (default: "mm10")
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_10X_h5(fname::String, dataset::String="/mm10"; unique_features::Bool=true)
    X, features, barcodes = _read_10X_h5(fname, dataset)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    convert_counts(X::AbstractMatrix, features::AbstractVector, barcodes::AbstractVector; unique_features::Bool=true)

Convert a count matrix and labels into its labeled representation

**Arguments**:

- `X`: a count matrix (features x barcodes)
- `features`: list of feature names
- `barcodes`: list of barcodes
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function convert_counts(X::AbstractMatrix, features::AbstractVector, barcodes::AbstractVector; unique_features::Bool=true)
    if !(eltype(X) <: Integer)
        @warn "count matrices should be integers, trying to convert from $(eltype(X))"
        X = convert(AbstractMatrix{Int64}, X)
    end

    if unique_features
        make_unique!(features, features)
    end

    NamedArray(copy(X'), (barcodes, features), (:cells, :features))
end

"""
    convert_counts(X::AbstractMatrix)

Convert a count matrix into its labeled representation by generating unique labels

**Arguments**:

- `X`: a count matrix (features x barcodes)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function convert_counts(X::AbstractMatrix)
    genes = [string("gene-", i) for i in 1:size(X,1)]
    barcodes = [string("cell-", i) for i in 1:size(X,2)]
    convert_counts(X, genes, barcodes; unique_features=false)
end

"""
    filter_counts(A::CountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)

Filter a count matrix, removing cells and features for which the metrics fall below the given threshold

**Arguments**:

    - `A`: the count matrix
    - `min_cells`: include features detected in at least this many cells
    - `min_features`: include cells where at least this many features are detected
    - `min_features_count`: threshold on the count for which a feature is marked "detected"
    - `min_umi`: include cells where the total of umi counts is at least this value

**Return value**:

The filtered matrix with cells and features removed
"""
function filter_counts(A::CountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)
    features_per_cell = vec(sum(A .> min_feature_count, dims=2))
    CI = (features_per_cell .>= min_features)

    if min_umi > 0
        CI .&= (vec(sum(A, dims=2)) .> min_umi)
    end

    if !all(CI)
        A = A[CI,:]
    end

    cells_per_feature = vec(sum(A .> 0, dims=1))
    FI = (cells_per_feature .>= min_cells)

    A[:, FI], CI, FI
    #  A[CI, FI] # faster but slightly different
end

"""
    filter_counts(A::NamedCountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)

Filter a labeled count matrix, removing cells and features for which the metrics fall below the given threshold

**Arguments**:

    - `A`: the count matrix
    - `min_cells`: include features detected in at least this many cells
    - `min_features`: include cells where at least this many features are detected
    - `min_features_count`: threshold on the count for which a feature is marked "detected"
    - `min_umi`: include cells where the total of umi counts is at least this value

**Return value**:

The filtered, labeled matrix with cells and features removed
"""
function filter_counts(A::NamedCountMatrix; min_cells=0, min_features=0, min_feature_count=0, min_umi=0)
    counts, CI, FI = filter_counts(A.array; min_cells=min_cells, min_features=min_features, min_feature_count=min_feature_count, min_umi=min_umi)
    barcodes, features = names(A)
    NamedArray(counts, (barcodes[CI], features[FI]), (:cells, :features))
end
