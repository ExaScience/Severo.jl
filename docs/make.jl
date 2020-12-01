using Documenter, Cell

makedocs(
    modules = [Cell],
    doctest = false,
    source = "src",
    build   = "build",
    clean = true,
    sitename = "Cell.jl",
    authors = "Tom Haber and contributors",
    pages = [
        "Home" => "index.md",
        "Showcases" => [
            "pbmc.md",
            "performance.md",
        ],
        "Library" => [
            "Public" => "public.md",
            "Internals" => "internals.md",
        ],
        "contributing.md",
    ]
)

deploydocs(;
    repo="github.imec.be/haber63/singlecell.git"
)
