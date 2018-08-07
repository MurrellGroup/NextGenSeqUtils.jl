using Documenter, NextGenSeqUtils
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/MurrellGroup/NextGenSeqUtils.jl.git",
    julia = "nightly",
    osx = "osx"
 
makedocs()
