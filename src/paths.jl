"""
Datatype for holding usearch and mafft file paths
"""
type Paths
    usearch::String
    mafft::String
end

"""
Default file paths.
"""
const PATHS = Paths("usearch", "mafft")
