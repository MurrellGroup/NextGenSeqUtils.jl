module NextGenSeqUtils 
    include("include_all.jl")
    
	export

    # align.jl

    usearch_filter,
    usearch_trim_fastq_with_phreds,

    nw_align,
    banded_nw_align,
    triplet_nw_align,
    local_align,
    kmer_seeded_align,
    triplet_kmer_seeded_align,
    loc_kmer_seeded_align,
    local_kmer_seeded_align,
    kmer_seeded_edit_dist,
    resolve_alignments,
    align_reference_frames,
    local_edit_dist,

    # hmm.jl

    viterbi_logs,
    trans_mat,
    obs_mat,
    initial_dist,
    get_obs_given_state,
    homopolymer_filter,
    markov_filter,
    forward_logs,
    backward_logs,
    forward_backward_logs,
    logsum,
    gen_seq_with_model,
    viterbiprint,

    # io.jl

    read_fasta,
    read_fasta_with_names,
    read_fasta_with_names_in_other_order,
    write_fasta,
    read_fastq,
    write_fastq,

    # kmers.jl

    KmerType,
    kmer_count,
    sparse_aa_kmer_count,
    corrected_kmer_dist,

    # orient.jl

    orient_strands,
    orient_to_refs,

    # paths.jl

    Paths,
    PATHS,

    # phreds.jl

    Phred,
    Prob,
    LogProb,
    MIN_PHRED,
    MAX_PHRED,
    phred_to_log_p,
    phred_to_p,
    p_to_phred,
    error_probs_to_phreds,
    quality_filter,
    length_vs_qual,
    qual_hist,

    # simulation.jl

    simple_gen_seq,
    simple_evolve,
    fixed_diff_evolve,
    pb_seq_sim,
    env_pb_seq_sim,

    # utils.jl

    print_fasta,
    degap,
    dash_count,
    single_gap,
    single_mod_three_gap,
    seq_details,
    print_rgb,
    reverse_complement,
    print_diffs,
    trim_ends_indices,
    translate_to_aa,
    generate_aa_seqs,
    filter_by_length,
    length_filter,
    concat_fastas,
    maxfreq,
    freq,
    sorted_freqs,
    freq_dict_print,
    nl43env

    # demux.jl
    demux_fastx
	
end # module
