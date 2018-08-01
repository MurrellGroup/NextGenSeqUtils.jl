# NextGenSeqUtils

## Synopsis

Contains a variety of functions for dealing with NGS data.

## Installation
```julia
Pkg.clone("NextGenSeqUtils")

```

## Set paths
```julia
using NextGenSeqUtils
```

## Files
**align.jl**
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


**hmm.jl**
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


**io.jl**
```julia
```

```julia
```

```julia
```

```julia
```

```julia
```

```julia
```

```julia
```

```julia
```

**kmers.jl**
>

**orient.jl**
>

**paths.jl**
>

**phreds.jl**
>

**simulation.jl**
>

**utils.jl**
>
