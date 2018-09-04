"""
    orient_strands(seqs::Array{String}, phreds::Union{Array{Vector{Phred},1},Nothing}, 
                   names::Union{Array{String,1},Nothing}, ref::String; k::Int=6)

Orients sequences (with phreds and names, which may be `nothing`) relative to a reference sequence.
`k` is kmer size for computing kmer vectors.
"""
function orient_strands(seqs::Array{String}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing}, ref::String; k::Int=6)
    kmers = [kmer_count(s, k) for s in seqs]
    ref_kmer = kmer_count(ref, k)
    rev_ref_kmer = kmer_count(reverse_complement(ref), k)
    #return kmer_orient_strands(seqs, ref, kmers, ref_kmer, rev_ref_kmer)

    # if this method crashes use the other return statement, it shouldnt work but it might
    return orient_strands(seqs, phreds, names, ref, kmers, ref_kmer, rev_ref_kmer)
end

"""
    orient_strands(seqs::Array{String}, ref::String; k::Int=6)

Orients sequences relative to a reference sequence.
`k` is kmer size for computing kmer vectors.
"""
orient_strands(seqs::Array{String}, ref::String; k::Int=6) = orient_strands(seqs, nothing, nothing, ref, k=k)[1]

"""
    orient_strands(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing}, 
                        ref::String, kmers::Array{Array{T,1},1}, ref_kmer::Array{T,1}, rev_ref_kmer::Array{T,1}) where {T <: Real}

Orients sequences with given kmer vectors relative to a reference sequence with given kmer vector.
"""
function orient_strands(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing}, 
                        ref::String, kmers::Array{Array{T,1},1}, ref_kmer::Array{T,1}, rev_ref_kmer::Array{T,1}) where {T <: Real}
    oriented = String[]
    newphreds = phreds==nothing ? nothing : Vector{Phred}[]
    newnames = names==nothing ? nothing : String[]

    for i in 1:length(seqs)
        vec = kmers[i]
        s = seqs[i]
        d_fwd = sqeuclidean(vec, ref_kmer)
        d_rev = sqeuclidean(vec, rev_ref_kmer)
        if d_fwd < d_rev
            push!(oriented, s)
            if phreds != nothing && length(phreds) > 0
                push!(newphreds, phreds[i])
            end
            if names != nothing && length(names) > 0
                push!(newnames, names[i])
            end
        else
            push!(oriented, reverse_complement(s))
            if phreds != nothing && length(phreds) > 0
                push!(newphreds, reverse(phreds[i]))
            end
            if names != nothing && length(names) > 0
                push!(newnames, names[i])
            end
        end
    end
    return oriented, newphreds, newnames
end

"""
    orient_to_refs_file(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, 
                        names::Union{Array{String,1},Nothing}, refspath::String)

Orients each sequence in `seqs` to nearest reference sequence in panel of references (`refspath`).
Distance determined by amino acid similarity (kmer vector dot prod).
"""
function orient_to_refs_file(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing}, refspath::String)
    if length(refspath) < 5
        return seqs
    end
    if refspath[end-5:end] == ".fasta"
        refs, _ = read_fasta_with_names_in_other_order(refspath, seqtype=String)
    elseif refspath[end-5:end] == ".fastq"
        refs, _, _ = read_fastq(refspath, seqtype=String)
    else
        error("Unsupported file type: $(infile[end-5:end])")
    end
    return orient_to_refs(seqs, phreds, names, refs)
end

"""
    orient_to_refs(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, 
                   names::Union{Array{String,1},Nothing}, refs::Array{String,1}; k::Int = 6)

Orients each sequence in `seqs` to nearest reference sequence in panel of references (`refs`).
Distance determined by amino acid similarity (kmer vector dot prod).
"""
function orient_to_refs(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing}, refs::Array{String,1}; k::Int = 6)
    if length(refs) == 0
        println("Not orienting: no reference sequences")
        return seqs
    end
    ref_kmers = [sparse_aa_kmer_count(rf, k) for rf in refs]
    rev_ref_kmers = [sparse_aa_kmer_count(reverse_complement(rf), k) for rf in refs]
    oriented = seqs[:]
    phredsexist = phreds != nothing && length(phreds) == length(seqs)
    namesexist = names != nothing && length(names) == length(seqs)
    function orientate(sq::String, ph::Vector{Phred}, name::String)
        best_match = 0
        best_match_i = 0
        rev = false
        vec = sparse_aa_kmer_count(sq, k)
        # compare to every reference, update min dists
        for i in 1:length(refs)
            match = dot(vec, ref_kmers[i])
            if match > best_match
                best_match = match
                best_match_i = i
                rev = false
            end
            match = dot(vec, rev_ref_kmers[i])
            if match > best_match
                best_match = match
                best_match_i = i
                rev = true
            end
        end
        if rev
            return reverse_complement(sq), reverse(ph), name
        elseif best_match_i > 0
            return sq, ph, name
        else
            # this probably happens because it does not share any kmers with any references
            # fix by shortening k?
            # make flag to throw these away?
            return sq, ph, name
        end            
    end
    oriented = pmap(i->orientate(oriented[i],
                                 phredsexist ? phreds[i] : Phred[],
                                 namesexist ? names[i] : ""),
                    [1:length(oriented)...])
    return [o[1] for o in oriented], phredsexist ? [o[2] for o in oriented] : nothing, namesexist ? [o[3] for o in oriented] : nothing
end

#=
"""Orient sequences without separate references. Coarse clusters and takes
refined consensus of largest cluster as reference"""
function orient_strands(seqs::Array{String, 1}; radius = 0.3, k = 4)
=#

