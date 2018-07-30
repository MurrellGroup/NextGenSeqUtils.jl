# Run `julia make.jl` in this folder to generate .html pages in a build/
# directory, then open index.md for docs

#push!(LOAD_PATH,"../src/")
using Documenter, NextGenSeqUtils

makedocs(
    format = :html,
    sitename = "NextGenSeqUtils.jl",
    modules = [NextGenSeqUtils],
    pages = [
        "index.md",
        "Documentation" => [
            "align.md",
            "hmm.md",
            "io.md",
            "kmers.md",
            "orient.md",
            "paths.md",
            "phreds.md",
            "simulation.md",
            "utils.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/MurrellGroup/NextGenSeqUtils.jl.git",
    julia = "0.6"
)

