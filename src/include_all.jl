using Reexport

# reexport gives the Levenshtein function in scope of file with "include" call
@reexport using Levenshtein
using BioSequences
import BioSequences: reverse_complement
using Distributions
using StatsBase
using Distances
using MultivariateStats

include("phreds.jl")
include("utils.jl")
include("io.jl")
include("kmers.jl")
include("simulation.jl")
include("align.jl")
include("hmm.jl")
include("orient.jl")
include("paths.jl")
include("wrappers.jl")
include("aliases.jl")
include("demux.jl")
