using Reexport

# reexport gives the Levenshtein function in scope of file with "include" call
#@reexport using Levenshtein
using BioSequences
import BioSequences: reverse_complement
using Distributions
using StatsBase
using Distances

#Uncomment when pull request to METADATA is completed by Multivar
#using MultivariateStats

using PyPlot
using LinearAlgebra

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
include("demux.jl")
include("FAD.jl")
include("evodist.jl")
include("levenshtein.jl")
