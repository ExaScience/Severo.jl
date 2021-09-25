# copyright imec - evaluation license - not for distribution

using Documenter, Severo

makedocs(
    modules = [Severo],
    doctest = false,
    source = "src",
    build   = "build",
    clean = true,
    sitename = "Severo.jl",
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
    repo="github.com/ExaScience/Severo.jl"
)
