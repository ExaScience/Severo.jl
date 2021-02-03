# copyright imec - evaluation license - not for distribution

import CSV
import DataFrames: DataFrame
import GZip
import HDF5
import HDF5: h5open, attributes, haskey, filename
import SparseArrays: nonzeros, rowvals, getcolptr

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
    ll = readline(mm, keep=false)
    while length(ll) == 0 || ll[1] == '%'
        ll = readline(mm, keep=false)
    end
    ll
end

function parseMM_comments(mm::IO)
    comments = []

    ll = readline(mm, keep=false)
    while length(ll) == 0 || ll[1] == '%'
        push!(comments, ll)
        ll = readline(mm, keep=false)
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
        ll = readline(mm, keep=false)

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

function _read_dge(fname::AbstractString; kw...)
    io = if endswith(fname, ".gz")
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end

    try
        X, comments = readMM(io; read_comments=true, kw...)

        genes = Vector{SubString}()
        sizehint!(genes, size(X,1))
        barcodes = Vector{SubString}()
        sizehint!(barcodes, size(X,2))

        for comment in comments
            if startswith(comment, "%%GENES")
                ss = split(comment, "\t")
                append!(genes, ss[2:end])
            elseif startswith(comment, "%%CELL_BARCODES")
                ss = split(comment, "\t")
                append!(barcodes, ss[2:end])
            end
        end

        copy(X'), genes, barcodes
    finally
        close(io)
    end
end

struct ParseError_10X <: Exception
    msg::AbstractString
end

function _read_10X_h5(fname::String, dataset::String="/mm10")
    h5open(fname, "r") do f
        feature_slot = if !haskey(attributes(f), "PYTABLES_FORMAT_VERSION")
            "/features/name"
        else
            "/gene_names"
        end

        if ! haskey(f, dataset)
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
            copy(X'), features, barcodes
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
        CSV.read(io, DataFrame; kw...)
    finally
        close(io)
    end
end

function _readlines(fname::AbstractString)
    io = if endswith(fname, ".gz")
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end

    try
        readlines(io)
    finally
        close(io)
    end
end

function _read_10X(dirname::AbstractString, gene_column::Int64=2)
end

function _read_10X(matrix_file::AbstractString, barcodes_file::AbstractString, feature_file::AbstractString, gene_column::Int64=2)
    if !(isfile(feature_file) && isfile(barcodes_file) && isfile(matrix_file))
        throw(ParseError_10X("Couldn't locate all components: $feature_file, $barcodes_file, $matrix_file"))
    end

    X = readMM(matrix_file)
    barcodes = _readlines(barcodes_file)
    features = readDelim(feature_file, header=false)[:,gene_column]
    copy(X'), features, barcodes
end

function _read_h5(fname::AbstractString, dataset::AbstractString="/counts")
    h5open(fname, "r") do f
        if ! haskey(f, dataset)
            throw(ArgumentError("Dataset $dataset does not exist in $fname"))
        end

        p = read(f, string(dataset, "/indptr"))
        i = read(f, string(dataset, "/indices"))
        x = read(f, string(dataset, "/data"))
        dim = read(f, string(dataset, "/shape"))
        features = read(f, string(dataset, "/rownames"))
        barcodes = read(f, string(dataset, "/colnames"))
        X = SparseMatrixCSC(dim[1], dim[2], p .+ 1, i .+ 1, x)
        copy(X'), features, barcodes
    end
end

struct ParseError_H5AD <: Exception
    msg::AbstractString
end

function read_h5ad_attr(attrs::HDF5.Attributes, desc::String, names::Vector{String})
    idx =  findfirst(x -> haskey(attrs, x), names)

    if idx === nothing
        throw(ArgumentError("Cannot read $desc information for count matrix in $(filename(attrs.parent))"))
    end

    read(attrs, names[idx])
end

function _read_h5ad(fname::AbstractString)
    h5open(fname, "r") do f
        if ! haskey(f, "X")
            throw(ArgumentError("Count data not found in $fname"))
        end

        a = attributes(f["X"])
        dim = read_h5ad_attr(a, "shape", ["shape", "h5sparse_shape"])
        format = read_h5ad_attr(a, "format", ["encoding-type", "h5sparse_format"])

        try
            p = read(f, "X/indptr") .+ 1
            i = read(f, "X/indices") .+ 1
            x = read(f, "X/data")

            X = if format == "csr" || format == "csr_matrix"
                X = SparseMatrixCSC(dim[2], dim[1], p, i, x)
                copy(X')
            else
                SparseMatrixCSC(dim[1], dim[2], p, i, x)
            end

            obs = read(f, "obs")
            barcodes = getindex.(obs, :index)

            var = read(f, "var")
            features = getindex.(var, :index)

            X, features, barcodes
        catch e
            if isa(e, ErrorException) # probably HDF5 error
                throw(ParseError_H5AD("Failed to load dataset $fname: $(e.msg)"))
            else
                rethrow(e)
            end
        end
    end
end

_keys(::Type{NamedTuple{names, types}}) where {names, types<:Tuple} = names

function _datatype(N::Type{<:NamedTuple})
    strtype = HDF5.HDF5Datatype(HDF5.h5t_copy(HDF5.H5T_C_S1))
    HDF5.h5t_set_cset(strtype, HDF5.H5T_CSET_UTF8)
    HDF5.h5t_set_size(strtype, HDF5.HDF5.H5T_VARIABLE)

    names = _keys(N)
    types = fieldtypes(N)

    size = 0
    for i in 1:nfields(types)
        T = types[i]
        data_type = if T == String
            strtype
        else
            HDF5.datatype(T)
        end
        size += sizeof(data_type)
    end

    dtype = HDF5.h5t_create(HDF5.H5T_COMPOUND, size)
    offset = 0
    for i in 1:nfields(types)
        T = types[i]
        data_type = if T == String
            strtype
        else
            HDF5.datatype(T)
        end

        HDF5.h5t_insert(dtype, String(names[i]), offset, data_type)
        offset += sizeof(data_type)
    end
    HDF5.HDF5Datatype(dtype)
end

function jl_to_hdf5(data::AbstractArray{<:NamedTuple}, i)
    N = eltype(data)
    T = fieldtype(N, i)

    if T == String
        ret = similar(data, Cstring)
        @inbounds for j in eachindex(data)
            ret[j] = Base.unsafe_convert(Cstring, data[j][i])
        end
        ret
    else
        ret = similar(data, T)
        @inbounds for j in eachindex(data)
            ret[j] = data[j][i]
        end
        ret
    end
end

function HDF5.write(parent::Union{HDF5.File, HDF5.Group}, name::String, data::AbstractArray{N}, plists::HDF5.Properties...) where {N<:NamedTuple}
    dtype = _datatype(N)
    dspace = HDF5.dataspace(data)

    strtype = HDF5.HDF5Datatype(HDF5.h5t_copy(HDF5.H5T_C_S1))
    HDF5.h5t_set_cset(strtype, HDF5.H5T_CSET_UTF8)
    HDF5.h5t_set_size(strtype, HDF5.HDF5.H5T_VARIABLE)

    try
        obj = HDF5.d_create(parent, name, dtype, dspace, plists...)

        try
            types = fieldtypes(N)
            names = _keys(N)
            for i in 1:nfields(types)
                T = types[i]
                data_type = if T == String
                    strtype
                else
                    HDF5.datatype(T)
                end
                tid = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(data_type))
                HDF5.h5t_insert(tid, String(names[i]), 0, data_type)
                HDF5.writearray(obj, tid, jl_to_hdf5(data, i))
                HDF5.h5t_close(tid)
            end
        finally
            close(obj)
        end
    finally
        close(dspace)
        close(dtype)
        close(strtype)
    end
end

function write_h5ad(fname::AbstractString, X::NamedCountMatrix)
    h5open(fname, "cw") do f
        try
            x = X.array
            write(f, "X/indptr", getcolptr(x) .- 1)
            write(f, "X/indices", rowvals(x) .- 1)
            write(f, "X/data", nonzeros(x))

            a = attributes(f["X"])
            a["encoding-type"] = "csc_matrix"
            a["shape"] = collect(size(x))

            nt = NamedTuple{(:index,),Tuple{String}}
            obs = map(n -> nt((n,)), names(X,1))
            write(f, "obs", obs)

            var = map(n -> nt((n,)), names(X,2))
            write(f, "var", var)
        catch e
            if isa(e, ErrorException) # probably HDF5 error
                throw(ParseError_H5AD("Failed to write to h5ad $fname: $(e.msg)"))
            else
                rethrow(e)
            end
        end

    end
end

function _read_csv(fname::AbstractString; unique_features::Bool=true)
    X = readDelim(fname)
    barcodes = names(X)[2:end]
    features = X[:,1]

    X = begin
        nz = 0
        for i in 2:size(X,2)
            nz += count(!iszero, X[!,i])
        end

        Tv = eltype(X[!,2])
        colptr = zeros(Int64, size(X, 2))
        rowval = Vector{Int64}(undef, nz)
        nzval = Vector{Tv}(undef, nz)
        colptr[1] = 1
        cnt = 1
        @inbounds for j in 2:size(X, 2)
            for i in 1:size(X, 1)
                v = X[i, j]
                if !iszero(v)
                    rowval[cnt] = i
                    nzval[cnt] = v
                    cnt += 1
                end
            end
            colptr[j] = cnt
        end
        SparseMatrixCSC(size(X, 1), size(X, 2) - 1, colptr, rowval, nzval)
    end

    copy(X'), features, barcodes
end

"""
    read_csv(dirname::AbstractString; unique_features=true)

Read count matrix from CSV

**Arguments**:

- `fname`: path to csv file
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_csv(dirname::AbstractString; unique_features::Bool=true)
    X, features, barcodes = _read_csv(dirname)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    read_csv(dirname::AbstractString; unique_features=true)

Read count matrix from digital gene expression (DGE) files

**Arguments**:

- `fname`: path to dge file
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_dge(dirname::AbstractString; unique_features::Bool=true)
    X, features, barcodes = _read_dge(dirname)
    convert_counts(X, features, barcodes, unique_features=unique_features)
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
function read_10X(dirname::AbstractString; gene_column::Int64=2, unique_features::Bool=true)
    if ! isdir(dirname)
        throw(ParseError_10X("Directory $dirname does not exist"))
    end

    is_v3 = isfile(joinpath(dirname, "features.tsv.gz"))
    feature_file = joinpath(dirname, is_v3 ? "features.tsv.gz" : "genes.tsv")
    barcodes_file = joinpath(dirname, is_v3 ? "barcodes.tsv.gz" : "barcodes.tsv")
    matrix_file = joinpath(dirname, is_v3 ? "matrix.mtx.gz" : "matrix.mtx")

    X, features, barcodes = _read_10X(matrix_file, barcodes_file, feature_file, gene_column)
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

function read_geo(prefix::AbstractString; gene_column::Int64=2, unique_features::Bool=true)
    is_gz = isfile(string(prefix, "matrix.mtx.gz"))
    postfix = is_gz ? ".gz" : ""

    feature_file = string(prefix, "genes.tsv", postfix)
    barcodes_file = string(prefix, "barcodes.tsv", postfix)
    matrix_file = string(prefix, "matrix.mtx", postfix)

    X, features, barcodes = _read_10X(matrix_file, barcodes_file, feature_file, gene_column)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    read_h5(fname::String, dataset::String="/mm10"; unique_features=true)

Read count matrix from hdf5 file.

**Arguments**:

- `fname`: path to hdf5 file
- `dataset`: name of dataset to load (default: "counts")
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_h5(fname::AbstractString, dataset::AbstractString="/counts"; unique_features::Bool=true)
    X, features, barcodes = _read_h5(fname, dataset)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    read_h5ad(fname::AbstractString, dataset::String="/mm10"; unique_features=true)

Read count matrix from hdf5 file as created by AnnData.py.
https://anndata.readthedocs.io/en/latest/fileformat-prose.html

**Arguments**:

- `fname`: path to hdf5 file
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_h5ad(fname::AbstractString; unique_features::Bool=true)
    X, features, barcodes = _read_h5ad(fname)
    convert_counts(X, features, barcodes, unique_features=unique_features)
end

"""
    read_data(path::AbstractString; unique_features=true)

Tries to identify and read a count matrix in any of the supported formats

**Arguments**:

- `fname`: path
- `unique_features`: should feature names be made unique (default: true)

**Returns values**:

Returns labeled sparse matrix containing the counts
"""
function read_data(path::AbstractString; unique_features::Bool=true)
    if !ispath(path)
        error("path $path does not point to anything")
    end

    f = if isdir(path)
        read_10X
    else
        if endswith(path, ".h5ad")
            read_h5ad
        elseif endswith(path, ".h5")
            read_10X_h5
        elseif endswith(path, ".csv")
            read_csv
        elseif endswith(path, ".dge.txt.gz") || endswith(path, ".dge.txt")
            read_dge
        else
            error("unknown format $path")
        end
    end

    f(path; unique_features=unique_features)
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
        X = convert(AbstractMatrix{Int32}, X)
    end

    if unique_features
        make_unique!(features, features)
    end

    NamedArray(X, (barcodes, features), (:cells, :features))
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
    barcodes = [string("cell-", i) for i in 1:size(X,1)]
    genes = [string("gene-", i) for i in 1:size(X,2)]
    convert_counts(X, genes, barcodes; unique_features=false)
end

