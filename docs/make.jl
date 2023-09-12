import Pkg
Pkg.add("Documenter")
using Documenter
using MolecularMinimumDistances
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [MolecularMinimumDistances],
    sitename = "MolecularMinimumDistances.jl",
    pages = [
        "Home" => "index.md",
        "Basic use" => "basic.md",
        "Advanced use" => "advanced.md",
        "Reference" => "reference.md",
        "Citation" => "citation.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/MolecularMinimumDistances.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
