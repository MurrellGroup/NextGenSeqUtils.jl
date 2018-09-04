using Reexport

# reexport gives the Levenshtein function in scope of file with "include" call
@reexport using Levenshtein
using BioSequences
import BioSequences: reverse_complement
using Distributions
using StatsBase
using Distances
#using MultivariateStats
#using NextGenSeqUtils
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
