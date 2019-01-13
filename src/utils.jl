#------Distance Matrix Functions--------
"""
    dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)

distances[i, j] is distance from distr1[i] to distr2[j].
"""
function dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)
    if distr1 == distr2
        return symmetric_dist_matrix(distr1, dist_met = dist_met)
    end
        
    distances = zeros(length(distr1), length(distr2))
    for i in 1:length(distr1)
        for j in 1:length(distr2)
            distances[i, j] = dist_met(distr1[i], distr2[j])
        end
    end
    return distances
end

"""
    symmetric_dist_matrix(distr; dist_met = kmer_seeded_edit_dist)

Similar to dist_matrix but where distr1 == distr2. 
Automatically called by dist_matrix if the two distributions are identical.
"""
function symmetric_dist_matrix(distr; dist_met = kmer_seeded_edit_dist)
    distances = zeros(length(distr), length(distr))
    for i in 1:length(distr)
        for j in i+1:length(distr)
            distances[i, j] = distances[j, i]  = dist_met(distr[i], distr[j])
        end
    end
    return distances
end


#-------OTHER UTILITY FUNCTIONS----------
"""
    print_fasta(seqs, names)

Prints fasta format to the terminal, for copypasting into alignment/blast etc.
"""
function print_fasta(seqs, names)
    for i in 1:length(seqs)
        println(">", names[i])
        println(seqs[i])
    end
end

"""
    degap(s::String)

Returns given string without '-' gap symbols.
"""
function degap(s::String)
    return replace(s, "-", "")
end

"""
    degap(s::DNASequence)

Returns given string without '-' gap symbols.
"""
function degap(s::DNASequence)
    return DNASequence(degap(String(s)))
end

"""
    dash_count(inStr::String)

Counts number of gap symbols '-' in given string.
"""
function dash_count(inStr::String)
    sum([char == '-' for char in inStr])
end

"""
    single_gap(str::String, ind::Int)

True if `str` has a single gap '-' at index `ind`, else false.
"""
function single_gap(str::String, ind::Int)
    return str[ind] == '-' && (ind == 1 || str[ind-1] != '-') && (ind == length(str) || str[ind+1] != '-')
end

"""
    single_mod_three_gap(str::String, ind::Int)

True if `str` has a gap length of 1 mod 3 at given index.
"""
function single_mod_three_gap(str::String, ind::Int)
    hi = ind
    while hi <= length(str) && str[hi] == '-'
        hi += 1
    end
    lo = ind
    while lo >= 1 && str[lo] == '-'
        lo -= 1
    end
    return (hi - lo - 1)%3 == 1
end

"""
    triple_gap(str::String, ind::Int)

True if index is within a gap of length a multiple of three, else false.
"""
function triple_gap(str::String, ind::Int)
    hi = ind
    while hi <= length(str) && str[hi] == '-'
        hi += 1
    end
    lo = ind
    while lo >= 1 && str[lo] == '-'
        lo -= 1
    end
    return (hi - lo - 1)%3 == 0
end

"""
    seq_details(fasta_path)

Gives names, sequences, error rates, and lengths from given filepath, which may end in '.fasta' or '.fastq'.
"""
function seq_details(fasta_path)
	if last(fasta_path) == 'a'
		names, seqs = read_fasta_with_names(fasta_path)
    	error_rates = [parse(Float64, split(i,['|', '='])[3]) for i in names]
	else
		seqs, errs, names = read_fastq(fasta_path)
		error_rates = [mean(phred_to_p(e)) for e in errs]
    end
    string_lengths = [length(i) for i in seqs]
    return names, seqs, error_rates, string_lengths
end

"""
    print_rgb(r, g, b, t)

Prints in colors `r`,`g`,`b` to terminal.
"""
function print_rgb(r, g, b, t)
    print("\e[1m\e[38;2;$r;$g;$b;249m", t)
end

"""
    reverse_complement(dna_string::String)

Returns the complement of the reverse of given nucleotide sequence.
"""
function reverse_complement(dna_string::String)
    return String(reverse_complement(DNASequence(dna_string)))
end

"""
    print_diffs(s1, s2; width=5, prefix="")

Prints two already aligned sequences with differences in color to terminal.
"""
function print_diffs(s1, s2; width=5, prefix="")
    if length(s1) != length(s2)
        error("Aligned strings are meant to be the same length.")
    end
    for i in 1:length(s1)
        if s1[i] != s2[i]
            println(prefix*"Position: ", i)
            # println(s1[max(i-width, 1):min(i+width, length(s1))])
            print_rgb(0, 0, 0, prefix*s1[max(i-width, 1):max(i-1, 1)])
            print_rgb(117, 41, 48, s1[i])
            print_rgb(0, 0, 0, s1[min(i+1, length(s1)):min(i+width, length(s1))])
            println("")

            print_rgb(0, 0, 0, prefix*s2[max(i-width, 1):max(i-1, 1)])
            print_rgb(117, 41, 48, s2[i])
            print_rgb(0, 0, 0, s2[min(i+1, length(s1)):min(i+width, length(s1))])
            println("")
        end
    end
end

"""
    trim_ends_indices(seq, ref; edge_reduction=0.1)

Align `seq` to `ref` with default low penalties for gaps on ends, 
and trim insertions on the ends of `seq`. 
Returns (start, stop) indices.
"""
function trim_ends_indices(seq, ref; edge_reduction=0.1)
    saln, raln = nw_align(seq, ref; edge_reduction=edge_reduction)
    x = 1
    while raln[x] == '-'
        x += 1
    end
    y = length(raln)
    while raln[y] == '-'
        y -= 1
    end
    n_end = length(raln) - y
    y = length(seq) - n_end
    return x, y
end

#----Amino Acids---

"""
    translate_to_aa(s::String)

Return amino acid string translation of nucleotide sequence using BioSequences conversion.
"""
function translate_to_aa(s::String)
    rna = convert(RNASequence, DNASequence(s))
    return string(translate(rna))
end

"""
    generate_aa_seqs(str::String)

Return sequence translated to amino acids in each reading frame
(returns three amino acid sequences)."""
function generate_aa_seqs(str::String)
    aa1 = translate_to_aa(str[1:length(str) - length(str)%3])
    aa2 = translate_to_aa(str[2:(length(str) - (length(str) - 1) % 3)])
    aa3 = translate_to_aa(str[3:(length(str) - (length(str) - 2) % 3)])
    return aa1, aa2, aa3
end

#----Length Filtering----

"""
    length_filter_inds(seqs::Array{String, 1}, minlength::Int, maxlength::Int)

Return indices of sequences in given array with length within given range (inclusive).
"""
function length_filter_inds(seqs::Array{String, 1}, minlength::Int, maxlength::Int)
    inds = map(s->(minlength <= length(s) <= maxlength), seqs)
    return inds
end

"""
    length_filter(seqs::Array{String, 1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing},
                  minlength::Int, maxlength::Int)

Filter sequences and corresponding names and phreds (which may be `nothing`) by length.
"""
function length_filter(seqs::Array{String, 1}, phreds::Union{Array{Vector{Phred},1},Nothing}, names::Union{Array{String,1},Nothing},
                           minlength::Int, maxlength::Int)
    inds = length_filter_inds(seqs, minlength, maxlength)
    seqs = seqs[inds]
    phreds = (phreds != nothing && length(phreds) > 0) ? phreds[inds] : phreds
    names = (names != nothing && length(names) > 0) ? names[inds] : names
    return seqs, phreds, names
end

"""
    length_filter(seqs::Array{String, 1}, minlength::Int64, maxlength::Int64)

Filter sequences by length.
"""
function length_filter(seqs::Array{String, 1}, minlength::Int64, maxlength::Int64)
    inds = length_filter_inds(seqs, minlength, maxlength)
    return seqs[inds]
end

"""
    filter_by_length(args...)

Deprecated. See `length_filter`.
"""
function filter_by_length(args...)
    Base.depwarn("Use length_filter instead.", :filter_by_length)
    return length_filter(args...)
end

"""
    concat_fastas(filepaths::Array{String, 1}, outfile::String)

Write contents of all given files to a single .fasta file.
"""
function concat_fastas(filepaths::Array{String, 1}, outfile::String)
    if length(outfile) < 7 || outfile[end-5:end] != ".fasta"
        error("Outfile must end in .fasta")
    end
    all_seqs, all_names = String[], String[]
    for (j, infile) in enumerate(infiles_list)
        fl = infile
        sqs, nms = read_fasta_with_names_in_other_order(fl, seqtype=String)
        for i in 1:length(sqs)
            push!(all_seqs, sqs[i])
            push!(all_names, nms[i])
        end
    end
    write_fasta(outfile, all_seqs, names=all_names)
end

#------Frequency Fxns------

"""
    maxfreq(vec)

Return the frequency of the most common element in `vec`.
"""
function maxfreq(vec)
    mo = mode(vec)
    proportionmap(vec)[mo]
end

"""
    freq(vec, elem)

Return the frequency of given element in given array; 
if the element is not present, return 0.0.
"""
function freq(vec, elem)
    map = proportionmap(vec)
    if haskey(map, elem)
        return map[elem]
        else return 0.0
    end
end

"""
    sorted_freqs(vec)

Return tuples of (freq, elem) of unique elements of `vec` 
in order of decreasing frequency.
"""
function sorted_freqs(vec)
    propDict = proportionmap(vec)
    seqkeys = keys(propDict)
    return reverse(sort([(propDict[k],k) for k in seqkeys]))
end

"""
    freq_dict_print(dictin; thresh=0)

Prints frequency:element of elements of `dictin` above given threshold, 
where `dictin` is a proportionmap of elements (see `proportionmap` in StatsBase).
"""
function freq_dict_print(dictin; thresh=0)
    for i in keys(dictin)
        if dictin[i]>thresh
            println(dictin[i], ":", i)
        end
    end
end

"""
    logsum(lga, lgb)

Compute numerically stable logsum. Returns `-Inf` if either input is `-Inf`.
"""
function logsum(lga, lgb)
    if lga == lgb == -Inf
        return -Inf
    end
    # log(x) + log(1 + exp(log(y) - log(x))) for x > y
    return lga > lgb ? (lga + log(1 + exp(lgb - lga))) : (lgb + log(1 + exp(lga - lgb)))
end


"""
	get_group_inds(arr, func)

Apply a function to a given set of objects, group by equal outputs, and return indices of elements in each group
"""
function get_group_inds(arr,func)
   groupings = [func(i) for i in arr]
   groups = union(groupings)
   groupinds = [find(groupings.==i) for i in groups]
   return groupinds
end


const nl43env = "GAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAGTATCAGCACTTGTGGAGATGGGGGTGGAAA"*
"TGGGGCACCATGCTCCTTGGGATATTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGT"*
"ACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATG"*
"TTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTT"*
"AACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATG"*
"TGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCG"*
"GGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGATAAGGTGCAG"*
"AAAGAATATGCATTCTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACAC"*
"CTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTG"*
"CGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACAT"*
"GGAATCAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGATGTAGTAATTAGATCTGC"*
"CAATTTCACAGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACA"*
"ACAATACAAGAAAAAGTATCCGTATCCAGAGGGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATG"*
"AGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATGCCACTTTAAAACAGATAGCTAGCAAATTAAGAGAACA"*
"ATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATT"*
"GTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAA"*
"GGGTCAAATAACACTGAAGGAAGTGACACAATCACACTCCCATGCAGAATAAAACAATTTATAAACATGTGGCAGGA"*
"AGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACTGGGCTGCTATTAA"*
"CAAGAGATGGTGGTAATAACAACAATGGGTCCGAGATCTTCAGACCTGGAGGAGGCGATATGAGGGACAATTGGAGA"*
"AGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGT"*
"GCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCG"*
"CAGCGTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGATATAGTGCAGCAGCAGAACAATTTGCTGAGG"*
"GCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAACAGCTCCAGGCAAGAATCCTGGCTGT"*
"GGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGC"*
"CTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATAACATGACCTGGATGGAGTGGGACAGAGAA"*
"ATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATT"*
"ATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCA"*
"TAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGA"*
"TATTCACCATTATCGTTTCAGACCCACCTCCCAATCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGG"*
"TGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTAGCACTTATCTGGGACGATCTGCGGAGCC"*
"TGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGG"*
"GGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAACTT"*
"GCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTATTACAAGCAGCTTATAGAGCTA"*
"TTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGG";


