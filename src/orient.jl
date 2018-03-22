"""Orients strands relative to a reference sequence."""
function orient_strands(seqs::Array{String}, ref::String; k::Int=6)
    kmers = [kmer_count(s, k) for s in seqs]
    ref_kmer = kmer_count(ref, k)
    rev_ref_kmer = kmer_count(reverse_complement(ref), k)
    return kmer_orient_strands(seqs, ref, kmers, ref_kmer, rev_ref_kmer)
end

"""Orients strands relative to a reference sequence."""
function orient_strands(seqs::Array{String}, ref::String, kmers::Array{T,1}, ref_kmer::Array{T,1}, 
                        rev_ref_kmer::Array{T,1}) where {T <: Real}
    oriented = String[]
    for i in 1:length(seqs)
        vec = kmers[i]
        s = seqs[i]
        d_fwd = sqeuclidean(vec, ref_kmer)
        d_rev = sqeuclidean(vec, rev_ref_kmer)
        if d_fwd < d_rev
            push!(oriented, s)
        else
            push!(oriented, reverse_complement(s))
        end
    end
    return oriented
end

"""Orients each sequence in `seqs` to nearest reference sequence in `refspath`.
Distance determined by amino acid similarity (kmer vector dot prod)"""
function orient_to_refs(seqs::Array{String, 1}, refspath::String)
    if length(refspath) < 5
        return seqs
    end
    if refspath[end-5:end] == ".fasta"
        refs = read_fasta_with_names_in_other_order(refspath, seqtype=String)
    elseif refspath[end-5:end] == ".fastq"
        refs = read_fastq(refspath, seqtype=String)
    else
        error("Unsupported file type: $(infile[end-5:end])")
    end
    return orient_to_refs(seqs, refs[1])
end

"""Orients each sequence in `seqs` to nearest reference sequence in `refs`.
Distance determined by amino acid similarity (kmer vector dot prod)"""
function orient_to_refs(seqs::Array{String, 1}, refs::Array{String, 1}; kmersize = 6)
    if length(refs) == 0
        println("Not orienting: no reference sequences")
        return seqs
    end
    ref_kmers = [sparse_aa_kmer_count(rf, kmersize) for rf in refs]
    rev_ref_kmers = [sparse_aa_kmer_count(reverse_complement(rf), kmersize) for rf in refs]
    oriented = seqs[:]
    function orientate(sq::String)
        best_match = 0
        best_match_i = 0
        rev = false
        vec = sparse_aa_kmer_count(sq, kmersize)
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
            return reverse_complement(sq)
        elseif best_match_i > 0
            return sq
        else
            # this probably happens because it does not share any kmers with any references
            # fix by shortening kmersize?
            # make flag to throw these away?
            return sq
        end            
    end
    oriented = pmap(i->orientate(oriented[i]), [1:length(oriented)...])
    return oriented
end

"""Orient sequences without separate references. Coarse clusters and takes
refined consensus of largest cluster as reference"""
#function orient_strands(seqs::Array{String, 1}; radius = 0.3, k = 4)
