# NextGenSeqUtils

## Synopsis

Contains a variety of functions for dealing with NGS data.

## Installation
```julia
Pkg.clone("https://github.com/MurrellGroup/NextGenSeqUtils.jl.git")

```

## Set paths
```julia
using NextGenSeqUtils
```

<a id='Files-1'></a>
# Files
**`FAD.jl`** &mdash; *File*
```julia
    FAD(seqs; alpha = 0.01, neigh_thresh = 1.0,method = 2,err_rate=0.02, phreds = nothing)

Exploits abundance and neighborhood information to denoise sequences without clustering on consensus calls. Primarily designed for short amplicons, high quality sequences, and good read-per-template coverage
```

**`align.jl`** &mdash; *File*
```julia
    nw_align(s1::String, s2::String; edge_reduction = 0.99)

Returns aligned strings using the Needleman-Wunch Algorithm (quadratic),
with end gaps penalized slightly less. edge_reduction is a multiplier (usually
less than one) on gaps on end of strings.
```

```julia
    nw_align(s1::String, s2::String, banded::Float64)

Wrapper for `nw_align` and `banded_nw_align`. A larger `banded` value makes alignment slower
but more accurate.
```

```julia
    add_to_band!(band, val, i::Int, j::Int, bandwidth::Int, dim_diff::Int)

Sets value in band at `i`, `j` to `val`, where `i` and `j` are in square matrix
coords, `dim_diff` = ncols - nrows
```

```julia
    get_band_val(band, i::Int, j::Int, bandwidth::Int, dim_diff::Int)

Returns value from band where `i` and `j` are in square matrix coords,
`dim_diff` = ncols - nrows
```

```julia
    in_band(i, j, bandwidth, dim_diff)

Checks if given square matrix coords are within width of band (ignore length)
```

```julia
    banded_nw_align(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)

Like nw_align, but sub quadratic by only computing values within a band around the center diagonal.
One 'band' of radius 3 = (4,1), (3,1), (2,1), (1,1), (1,2), (1,3), (1,4), aka upside-down L shape.
band_coeff = 1 is sufficient to get same alignments as nw_align for 10% diverged sequences ~97% of the time;
increase this value for more conservative alignment with longer computation time.
Radius of band = `bandwidth` = `band_coeff` * sqrt(avg seq length)
```

```julia
    triplet_nw_align(s1::String, s2::String; edge_reduction = 0.99, boundary_mult = 2)

Returns alignment of two sequences where `s1` is a reference with reading frame to be preserved and `s2` is a query sequence.
`boundary_mult` adjusts penalties for gaps preserving the reading frame of `s1`.
This usually works best on range 0 to 3, higher values for more strongly enforced gaps aligned on
reference frame (divisible-by-3 indices)
```

```julia
    local_align(ref::String, query::String; mismatch_score = -1,
            match_score = 1, gap_penalty = -1,
            rightaligned=true, refend = false)

Aligns a query sequence locally to a reference. If true, `rightaligned` keeps the
right ends of each sequence in final alignment- otherwise they are trimmed;
`refend` keeps the beginning/left end of `ref`.
If you want to keep both ends of both strings, use nw_align.
For best alignments use the default score values.
```

```julia
    unique_key(dicto::Dict{String, Int}, keyo::String, indo::Int)

    Sets value of key in dictionary to given index, unless the key already exists,
    in which case value is set to -1.
```

```julia
    sorted_matches(s1, s2, wordlength, skip, aligncodons)

    Returns a list of indices of matches of unique words between `s1` and `s2`,
sorted by index of `s1`.
```

```julia
    sorted_aa_matches(str1, str2, wordlength)

    Returns a list of indices of matches of unique words between `s1` and `s2`,
    sorted by index of `s1`. Words match by amino acid encoding, in any reading frame.    
```

```julia
    matches_are_inconsistent(matches)

    Returns true if indices of word matches in second sequence (second column of
    `matches`) is strictly increasing, else returns false.

```

```julia
    clean_matches(matches, wordlength, skip)

Some edge cases commonly arise where a kmer match starts before
the previous kmer match ends, but the two sequences still
mismatch.

To handle these, we make sure that either the word match in both sequences is spaced
`skip` after the previous word match, or both matches have spacing greater than
`wordlength`.
```


```julia
    merge_overlapping(matches, wordlength, skip)

Returns ranges of overlapping word matches.
```

```julia
    get_matches(s1, s2, clean, range_inds, wordlength)

Returns actual word matches from indices of the matches.
```

```julia
    get_mismatches(s1, s2, matches, range_inds, wordlength)
Returns actual word mismatches (intervals between word matches) from indices of matches.
```

```julia
    longest_incr_subseq(arr::Array{Int, 2})

Find longest increasing subsequence of second column of given array.
Used to resolve bad orders of word matches while preserving as many matches as possible.
Because the first column (match locations in the first sequence) is sorted, corresponding
matches in the second column must be sorted, so we get the maximum number of such matches.
```

```julia
    kmer_seeded_align(s1::String, s2::String;
                  wordlength = 30,
                  skip = 10,
                  aligncodons = false,
                  banded = 1.0,
                  debug::Bool = false)

Returns aligned strings, where alignment is first done with larger word matches and
then (possibly banded) Needleman-Wunsch on intermediate intervals.
`skip` gives a necessary gap between searched-for words in `s1`.
For best results, use the default `wordlength` and `skip` values.
See `nw_align` for explanation of `banded`.
```


```julia
triplet_kmer_seeded_align(s1::String, s2::String;
                          wordlength = 30,
                          skip = 9,
                          boundary_mult = 2,
                          alignedcodons = true,
                          debug::Bool=false)

Returns aligned strings, where alignment is first done with word matches and
then Needleman-Wunsch on intermediate intervals, prefering to preserve the
reading frame of the first arg `s1`.
`skip` gives a necessary gap between searched-for words in `s1`.
For best results, use the default `wordlength` and `skip` values.
See `triplet_nw_align` for explanation of `boundary_mult`.
```

```julia
function local_kmer_seeded_align(s1::String, s2::String;
                                 wordlength = 30,
                                 skip = 10,
                                 trimpadding = 100,
                                 debug::Bool=false)

Returns locally aligned strings, where alignment is first done with word matches and
then Needleman-Wunsch on intermediate intervals.

`s1` is a reference to align to, and `s2` is a query to extract a local match from.
`s2` may be trimmed or expanded with gaps.
Before locally aligning ends of sequences, the ends of `s2` are trimmed to length
`trimpadding` for faster alignment. Increasing this will possible increase alignment accuracy
but effect runtime.
`skip` gives a necessary gap between searched-for words in `s1`.
For best results, use the default `wordlength` and `skip` values.
```

```julia
kmer_seeded_edit_dist(s1::String , s2::String;
                      wordlength = 30,
                      skip = 5,
                      aa_matches = false)

Computes levenshtein edit distance with speedups from only computing the dp scoring matrix between word matches.
If aa_matches = true, will attempt to find amino acid matches in any reference frame,
and add the nucleotide Hamming distance of these matches to Levenshtein distances of mismatches.
`skip` gives a necessary gap between searched-for words in `s1`.
For best results, use the default `wordlength` and `skip` values.
```

```julia
resolve_alignments(ref::String, query::String; mode = 1)

Called on aligned strings. Resolves `query` with respect to `ref`.
`mode` = 1 for resolving single indels, `mode` = 2 for resolving single indels and codon insertions in query.
```

```julia
align_reading_frames(clusters; k = 6, thresh = 0.03, verbose = false)

Takes `clusters` = [consensus_sequences, cluster_sizes], chooses references
out of consensuses that do not have stop codons in the middle, and makes all
consensus sequence reading frames agree. Returns resolved consensus seqs
(`goods`) along with filtered out consensus seqs that are >`thresh` divergent from
nearest reference (`bads`).
`k` = kmer size for computing kmer vectors of sequences.
```

```julia
local_edit_dist(s1::String, s2::String)

Returns the edit distance between two sequences after local alignment
```


**`hmm.jl`** &mdash; *File*
```julia
    viterbi_logs(observations_given_states::Array{Float64, 2}, transitions::Array{Float64, 2}, initials::Array{Float64})

Return viterbi path and log probability for that path. Takes logs of matrices.
`observations_given_states` has # rows = # states, # columns = # steps of markov process.
```

```julia
trans_mat(; uniform_cycle_prob = 0.9999999999, homopoly_cycle_prob = 0.98)

Return 5x5 transition matrix with given transition probabilities
```

```julia
obs_mat(; homopoly_prob = 0.99)

Return 5x4 observation matrix with given probabilities.
```

```julia
initial_dist(; uniform_state = 0.99)

Initial state distributions: see `trans_mat` for state descriptions. 5x1 vector.
```

```julia
get_obs_given_state(observation_matrix::Array{Float64,2}, observation_seq::String)

Populate 5xT matrix with likelihood of observation at each of T time steps given each
nucleotide in `observation_seq`.
```

```julia
homopolymer_filter(seqs::Array{String,1}, phreds, names;
                   transmat = nothing, obsmat = nothing,
                   initialdist = nothing)

Filter sequences with "bad" sections in the middle -- abnormally long runs of a single base, using
viterbi alg inference. If this homopolymer occurs on one end of a sequence, keeps sequence and trims
homopolymer region off end. `phreds` and/or `names` may be `nothing`, in which case `nothing` is returned
for the respective field. If transition, observation, initial distribution matrices not provided (default)
then the default values from the respective constructors are used.
```

```julia
homopolymer_filter(sourcepath::String, destpath::String;
                   transmat = nothing, obsmat = nothing,
                   initialdist = nothing, format="fastq")

Filter sequences with "bad" sections in the middle -- abnormally long runs of a single base, using
viterbi alg inference. If this homopolymer occurs on one end of a sequence, keeps sequence and trims
homopolymer region off end. Takes a file path for each of a source file of type `format` (which may be "fasta" or "fastq")
and a destination file is written which is the same file type.
```

```julia
forward_logs(observations_given_states::Array{Float64, 2}, transitions::Array{Float64, 2}, initials::Array{Float64})

Compute logs of forward scores. Takes logs of matrices as inputs.
`observations_given_states` has # rows = # states, # columns = # steps of markov process.
```

```julia
backward_logs(observations_given_states::Array{Float64, 2}, transitions::Array{Float64, 2}, initials::Array{Float64})

Compute logs of backward scores and individual posterior probabilities. Takes logs of matrices as inputs.
`observations_given_states` has # rows = # states, # columns = # steps of markov process.
```

```julia
forward_backward_logs(observations_given_states::Array{Float64, 2}, transitions::Array{Float64, 2}, initials::Array{Float64})

Compute logs of forward-backward scores and individual probabilities. Takes logs of matrices as inputs.
`observations_given_states` has # rows = # states, # columns = # steps of markov process.
```

```julia
gen_seq_with_model(n::Int, trans_mat, obs_mat, initial_dists)

Generate a sequence given a markov model. May create bad reads (long runs of a single base).
```

```julia
viterbiprint(s::String)

Draw flagged sites with capital letters, safe sites with lowercase.
```


**`io.jl`** &mdash; *File*
```julia
read_fasta_records(filename)

Read .fasta file contents.
```

```julia
read_fasta(filename; seqtype=String)

Read .fasta file contents, parse, and return sequences as type `seqtype`
```

```julia
read_fasta_with_names(filename; seqtype=String)

Read .fasta file contents, parse, and return names and sequences as type `seqtype`.
```

```julia
read_fasta_with_names_in_other_order(filename; seqtype=String)

Read .fasta file contents, parse, and return sequences as type `seqtype` with names.
Now that's what I call convenience.
```

```julia
write_fasta(filename::String, seqs; names = String[])

Write given `seqs` and optional `names` to a .fasta file with given filepath.
```

```julia
read_fastq_records(filename)

Read .fastq file contents.
```

```julia
read_fastq(filename; seqtype=String)

Read .fastq file contents, parse, and return sequences as `seqtype` type, phreds, and names.
```

```julia
write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                 names=String[], DNASeqType = false)

Write given sequences, phreds, names to .fastq file with given file path.
If `names` not provided, gives names 'seq_1', etc.
```

**`kmers.jl`** &mdash; *File*
```julia
kmer_count(str::String, k::Int)

Count kmers of size `k` in string, return array with bins from "A...A" to "T...T"
and value in each bin corresponding to number of occurences of that kmer within `str`.
```

```julia
sparse_aa_kmer_count(str::String, k::Int)

Counts amino acid kmers in string (in all reading frames).
k = length of kmer in amino acids.
This is sparse and kinda slow.
```

```julia
corrected_kmer_dist(kmers1::Array, kmers2::Array; k = nothing)

Compute distance between kmer vectors that is corrected towards edit distance for small differences.
The default for `k` is `k = Int(log(4, length(kmers1)))`
```

```julia
corrected_kmer_dist(k::Int)

Returns function that computes distance between `k`-mer vectors
that is corrected towards edit distance for small differences.
```

**`orient.jl`** &mdash; *File*
```julia
orient_strands(seqs::Array{String}, phreds::Union{Array{Vector{Phred},1},Void},
               names::Union{Array{String,1},Void}, ref::String; k::Int=6)

Orients sequences (with phreds and names, which may be `nothing`) relative to a reference sequence.
`k` is kmer size for computing kmer vectors.
```

```julia
orient_strands(seqs::Array{String}, ref::String; k::Int=6)

Orients sequences relative to a reference sequence.
`k` is kmer size for computing kmer vectors.
```

```julia
orient_strands(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Void}, names::Union{Array{String,1},Void},
                    ref::String, kmers::Array{Array{T,1},1}, ref_kmer::Array{T,1}, rev_ref_kmer::Array{T,1}) where {T <: Real}

Orients sequences with given kmer vectors relative to a reference sequence with given kmer vector.
```

```julia
orient_to_refs_file(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Void},
                    names::Union{Array{String,1},Void}, refspath::String)

Orients each sequence in `seqs` to nearest reference sequence in panel of references (`refspath`).
Distance determined by amino acid similarity (kmer vector dot prod).
```

```julia
orient_to_refs(seqs::Array{String,1}, phreds::Union{Array{Vector{Phred},1},Void},
               names::Union{Array{String,1},Void}, refs::Array{String,1}; k::Int = 6)

Orients each sequence in `seqs` to nearest reference sequence in panel of references (`refs`).
Distance determined by amino acid similarity (kmer vector dot prod).
```

**`phreds.jl`** &mdash; *File*
```julia
phred_to_log_p(x)

Conversion from phred value to log probability value.
```

```julia
phred_to_p(q::Phred)

Conversion from phred value to real probability value.
```

```julia
phred_to_p(x::Vector{Phred})

Conversion from phred values to real probability values in vector.
```

```julia
p_to_phred(p::Prob)

Conversion from real probability value to phred value.
```

```julia
p_to_phred(p::Prob)

Conversion from real probability values to phred values in vector.
```

```julia
error_probs_to_phreds(eps::Vector{Float64})

Conversion from error probability values to phred values in vector.
```

```julia
quality_filter(infile, outfile=join(split(infile, ".")[1:end-1], ".") * ".filt.fasta";
               errorRate = 0.01, minLength = 100, labelPrefix = "seq", errorOut = true)

Writes file with sequences from input file that have all sites within error rate margin.
```

```julia
quality_filter(seqs::Array{String, 1}, scores::Array, names::Array{String, 1}; errorRate=0.01, minLength=0)

Returns sub arrays with sequences that have site-wise error all within given margin.
```

```julia
quality_filter_inds(seqs::Array{String, 1}, scores::Array; errorRate=0.01,minLength=100)

Returns indices of sequences with site-wise error within given margin.
```

```julia
length_vs_qual(fasta_path; plot_title = "Length Vs Errors")

Creates plot of lengths of sequences vs. mean error rates of sequences.
```

```julia
qual_hist(fasta_path; plot_title = "Quality Hist")

Creates histogram of mean quality scores per sequence.
```


**`simulation.jl`** &mdash; *File*
```julia
flip(p, t, f)

Flip a biased coin.
```

```julia
simple_gen_seq(n::Int)

Generates a random uniform sequence of nucleotides.
```

```julia
simple_evolve(refseq, err_rate)

Evolves a sequence, uniformly.
```

```julia
fixed_diff_evolve(template::String, n_diffs::Int64)

Creates a fixed number of mutations of a sequence.
```

```julia
run_length_encode(x::String)

Run length encoding of a string.
```

```julia
run_length_decode(chars, lengths)

Run length decoding of a string.
```

```julia
pb_error_inflation(old_length::Int64)

Homopolymer length-to-error-rate scaling based on PacBio sequencing.

This function may need to be tweaked to approximate how the error
rate increases with HP length.
```

```julia
length_error_func(old_length::Int64; rate = 0.002)

Takes a true homopolymer length, and returns an observed homopolymer length based on PacBio sequencer error model.
```

```julia
pb_seq_sim(template::String, rate::Float64; with_qvs = false)

Performs a sequence simulation from a template, specifying a target error rate, based on PacBio sequencing error model.
```

```julia
env_error_rates(n)

Draws from the error rate distribution typically seen in P5 envelope
sequence data.
```

```julia
env_pb_seq_sim(template::String, n::Int64; with_qvs = false)

Simulated PacBio reads from amplicons that have envelope-like error profiles.
```

**`utils.jl`** &mdash; *File*
```julia
dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)

distances[i, j] is distance from distr1[i] to distr2[j].
```

```julia
symmetric_dist_matrix(distr; dist_met = kmer_seeded_edit_dist)

Similar to dist_matrix but where distr1 == distr2.
Automatically called by dist_matrix if the two distributions are identical.
```

```julia
print_fasta(seqs, names)

Prints fasta format to the terminal, for copypasting into alignment/blast etc.
```

```julia
degap(s::String)

Returns given string without '-' gap symbols.
```

```julia
degap(s::DNASequence)

Returns given string without '-' gap symbols.
```

```julia
dash_count(inStr::String)

Counts number of gap symbols '-' in given string.
```

```julia
single_gap(str::String, ind::Int)

True if `str` has a single gap '-' at index `ind`, else false.
```

```julia
single_mod_three_gap(str::String, ind::Int)

True if `str` has a gap length of 1 mod 3 at given index.
```

```julia
triple_gap(str::String, ind::Int)

True if index is within a gap of length a multiple of three, else false.
```

```julia
seq_details(fasta_path)

Gives names, sequences, error rates, and lengths from given filepath, which may end in '.fasta' or '.fastq'.
```

```julia
print_rgb(r, g, b, t)

Prints in colors `r`,`g`,`b` to terminal.
```

```julia
reverse_complement(dna_string::String)

Returns the complement of the reverse of given nucleotide sequence.
```

```julia
print_diffs(s1, s2; width=5, prefix="")

Prints two already aligned sequences with differences in color to terminal.
```

```julia
trim_ends_indices(seq, ref; edge_reduction=0.1)

Align `seq` to `ref` with default low penalties for gaps on ends,
and trim insertions on the ends of `seq`.
Returns (start, stop) indices.
```

```julia
translate_to_aa(s::String)

Return amino acid string translation of nucleotide sequence using BioSequences conversion.
```

```julia
generate_aa_seqs(str::String)

Return sequence translated to amino acids in each reading frame
(returns three amino acid sequences).
```

```julia
length_filter_inds(seqs::Array{String, 1}, minlength::Int, maxlength::Int)

Return indices of sequences in given array with length within given range (inclusive).
```

```julia
length_filter(seqs::Array{String, 1}, phreds::Union{Array{Vector{Phred},1},Void}, names::Union{Array{String,1},Void},
              minlength::Int, maxlength::Int)

Filter sequences and corresponding names and phreds (which may be `nothing`) by length.
```

```julia
length_filter(seqs::Array{String, 1}, minlength::Int64, maxlength::Int64)

Filter sequences by length.
```

```julia
filter_by_length(args...)

Deprecated. See `length_filter`.
```

```julia
concat_fastas(filepaths::Array{String, 1}, outfile::String)

Write contents of all given files to a single .fasta file.
```

```julia
maxfreq(vec)

Return the frequency of the most common element in `vec`.
```

```julia
freq(vec, elem)

Return the frequency of given element in given array;
if the element is not present, return 0.0.
```

```julia
sorted_freqs(vec)

Return tuples of (freq, elem) of unique elements of `vec`
in order of decreasing frequency.
```

```julia
freq_dict_print(dictin; thresh=0)

Prints frequency:element of elements of `dictin` above given threshold,
where `dictin` is a proportionmap of elements (see `proportionmap` in StatsBase).
```

```julia
logsum(lga, lgb)

Compute numerically stable logsum. Returns `-Inf` if either input is `-Inf`.
```

**`evodist.jl`** &mdash; *File*
```julia
    get_pis(sequences)
```

```julia
    get_transition_mat(seq1, seq2; include_dash_seqs = true)
```

```julia
    estimate_distance(seq1, seq2)
```

