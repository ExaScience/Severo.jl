# copyright imec - evaluation license - not for distribution

import HTTP
import GZip
import Tar

data_location(collection::AbstractString="") = joinpath(@__DIR__, "..", "data", collection)

abstract type DataSet end

struct FileDataSet <: DataSet
    name::String
    filename::String
    url::String
end

struct TarDataSet <: DataSet
    description::String
    name::String
    url::String
end

struct DataCollection
    short::String
    description::String
    doi::String
    link::String
    datasets::Vector{DataSet}
end

function preprocess(collections...)
    d = Dict{String, DataCollection}()
    for col in collections
        push!(d, col.short => col)
    end
    d
end

const datacollections = preprocess(
    DataCollection("Svensson2019", "Negative control data for drop-seq", "10.1101/582064",
        "https://figshare.com/projects/Zero_inflation_in_negative_control_data/61292", [
        FileDataSet("gemcode", "zheng_gemcode_control.h5ad", "https://ndownloader.figshare.com/files/14634407"),
        FileDataSet("chromium", "svensson_chromium_control.h5ad", "https://ndownloader.figshare.com/files/14634410"),
        FileDataSet("indrops", "klein_indrops_control_GSM1599501.h5ad", "https://ndownloader.figshare.com/files/14634416"),
        FileDataSet("dropseq", "macosko_dropseq_control.h5ad", "https://ndownloader.figshare.com/files/14634488")
    ]),
    DataCollection("Zheng2017", "Massively parallel digital transcriptional profiling of single cells", "10.1038/ncomms14049",
        "https://support.10xgenomics.com/single-cell-gene-expression/datasets", [
        TarDataSet("CD19+ B Cells", "b_cells", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD14+ Monocytes", "cd14_monocytes", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_raw_gene_bc_matrices.tar.gz"),
        TarDataSet("CD34+ Cells", "cd34", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd34/cd34_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD4+ Helper T Cells", "cd4_t_helper", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/cd4_t_helper_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD4+/CD25+ Regulatory T Cells", "regulatory_t", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/regulatory_t/regulatory_t_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD4+/CD45RA+/CD25- Naive T cells", "naive_t", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_t/naive_t_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD4+/CD45RO+ Memory T Cells", "memory_t", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/memory_t/memory_t_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD56+ Natural Killer Cells", "cd45_nk", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/cd56_nk_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD8+ Cytotoxic T cells", "cytotoxic_t", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cytotoxic_t/cytotoxic_t_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("CD8+/CD45RA+ Naive Cytotoxic T Cells", "naive_cytotoxic", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/naive_cytotoxic/naive_cytotoxic_filtered_gene_bc_matrices.tar.gz"),
        TarDataSet("Fresh 68k PBMCs (Donor A)", "pbmc_68k", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz")
    ]),
    DataCollection("1.3M_nn", "1.3 Million Brain Cells from E18 Mice", "",
        "https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons", [
            FileDataSet("full", "1M_neurons_full.h5", "https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5")
            FileDataSet("sampled", "1M_neurons_sampled.h5", "https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_neuron20k.h5")
        ]),
    DataCollection("Saunders2018", "Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain", "10.1016/j.cell.2018.07.028",
        "http://dropviz.org/", [
            FileDataSet("cerebellum", "cerebellum.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz"),
            FileDataSet("entopeduncular", "entopeduncular.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60EntoPeduncular.raw.dge.txt.gz"),
            FileDataSet("frontal_cortex", "frontal_cortex.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz"),
            FileDataSet("globus_pallidus", "globus_pallidus.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60GlobusPallidus.raw.dge.txt.gz"),
            FileDataSet("hippocampus", "hippocampus.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"),
            FileDataSet("posterior_cortex", "posterior_cortex.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz"),
            FileDataSet("stratium", "stratium.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Striatum.raw.dge.txt.gz"),
            FileDataSet("substantia_nigra", "substantia_nigra.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60SubstantiaNigra.raw.dge.txt.gz"),
            FileDataSet("thalamus", "thalamus.raw.dge.txt.gz", "https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Thalamus.raw.dge.txt.gz")
        ]),
    DataCollection("PBMC", "Peripheral blood mononuclear cells (PBMCs) from a healthy donor", "",
        "https://support.10xgenomics.com/single-cell-gene-expression/datasets", [
            TarDataSet("33k PBMCs from a Healthy Donor", "33k", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc33k/pbmc33k_filtered_gene_bc_matrices.tar.gz"),
            TarDataSet("3k PBMCs from a Healthy Donor", "3k", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"),
            TarDataSet("6k PBMCs from a Healthy Donor", "6k", "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc6k/pbmc6k_filtered_gene_bc_matrices.tar.gz")
        ])
)

list_collections() = collect(values(datacollections))

function list_datasets(collection::AbstractString)
    col = get(datacollections, collection, nothing)
    if col === nothing
        error("collection $collection not found")
    end

    col.datasets
end

show(io::IO, ds::FileDataSet) = print(io, "Dataset $(ds.name)")
show(io::IO, ds::TarDataSet) = print(io, "Dataset $(ds.name): $(ds.description)")

function show(io::IO, collection::DataCollection)
    println(io, "Collection $(collection.short)")
    println(io, "  $(collection.description)")
    println(io, "  DOI: $(collection.doi)")
    println(io, "  WEB: $(collection.link)")
    println(io)
    for ds in collection.datasets
        println(io, "    ", ds)
    end
end

function download(url::AbstractString, location::AbstractString)
    resp = HTTP.open("GET", url) do stream
        resp = HTTP.startread(stream)

        total_bytes = parse(Int64, HTTP.getkv(resp.headers, "Content-Length", 0))

        if resp.status != 200
            return
        end

        @info "Downloading $url" size=total_bytes location=location

        bytes_downloaded = open(location, "w") do fh
            bytes = 0
            while(!eof(stream))
                bytes += write(fh, readavailable(stream))
            end
            bytes
        end

        if total_bytes != 0 && bytes_downloaded != total_bytes
            rm(location)
            error("download of $url failed: length differs $total_bytes != $bytes_downloaded")
        end

        @info "Finished downloading" bytes=bytes_downloaded
    end

    if resp.status != 200
        error("downloading $url has failed: status=$(resp.status)")
    end
end

function download(url::AbstractString, io::IO)
    resp = HTTP.open("GET", url) do stream
        resp = HTTP.startread(stream)

        total_bytes = parse(Int64, HTTP.getkv(resp.headers, "Content-Length", 0))

        if resp.status != 200
            return
        end

        @info "Downloading $url" size=total_bytes

        bytes_downloaded = 0
        while(!eof(stream))
            bytes_downloaded += write(io, readavailable(stream))
        end

        if total_bytes != 0 && bytes_downloaded != total_bytes
            error("download of $url failed: length differs $total_bytes != $bytes_downloaded")
        end

        @info "Finished downloading" bytes=bytes_downloaded
    end

    if resp.status != 200
        error("downloading $url has failed: status=$(resp.status)")
    end
end

function maybe_download(ds::FileDataSet, path::AbstractString)
    location = joinpath(path, ds.filename)
    if ! isfile(location)
        download(ds.url, location)
    end

    location
end

function extract_tarball(tar::IO, location::AbstractString)
    mktempdir() do path
        list = Tar.Header[]
        Tar.extract(h -> (push!(list, h); true), tar, path)
        for hdr in list
            if hdr.type == :file
                mv(joinpath(path, hdr.path), joinpath(location, basename(hdr.path)))
            end
        end
    end
end

function maybe_download(ds::TarDataSet, path::AbstractString)
    location = joinpath(path, ds.name)
    if ! ispath(location)
        mkdir(location)

        mktemp() do path, io
            download(ds.url, io)
            seekstart(io)
            tar = GZip.gzdopen(io)
            try
                extract_tarball(tar, location)
            finally
                close(tar)
            end
        end
    end

    location
end

function dataset(collection::AbstractString, dataset::AbstractString)
    col = get(datacollections, collection, nothing)
    if col === nothing
        error("collection $collection not found")
    end

    idx = findfirst(ds -> ds.name == dataset, col.datasets)
    if idx === nothing
        error("dataset ($dataset) not found in collection ($collection)")
    end

    ds = col.datasets[idx]
    path = joinpath(data_location(), col.short)

    if ! ispath(path)
        mkdir(path)
    end

    read_data(maybe_download(ds, path))
end
