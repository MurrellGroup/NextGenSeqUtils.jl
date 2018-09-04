"""
Datatype for holding usearch and mafft file paths
"""
mutable struct Paths
    usearch::String
    mafft::String
end

"""
Default file paths.
"""
const PATHS = Paths("usearch", "mafft")
