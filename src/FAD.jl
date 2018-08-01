"""
    FAD(seqs; alpha = 0.01, neigh_thresh = 1.0,method = 2,err_rate=0.02, phreds = nothing)

Exploits abundance and neighborhood information to denoise sequences without clustering
on consensus calls. Primarily designed for short amplicons, high quality sequences,
and good read-per-template coverage
"""

function FAD(seqs; alpha = 0.01, neigh_thresh = 1.0,method = 2,err_rate=0.02, phreds = nothing)
    if !(method in [1,2,3])
        warn("Please pick a method in 1,2,3")
        return [],[]
    end
    #Method list:
    #1: No statistical considerations. All smaller relatives get eaten.
    #2: Use Poisson test to decide to add not-too-small relatives back in.
    #3: Use homopolymer test to skew Poisson rate for adding smaller relatives back in.

    #Need to mod this to work for .fasta files. Perhaps use the number of pairs to approximate the number of error free sequences.
    if phreds != nothing
        warn("Switching to method 1, phreds missing")
        method = 1
        #seqs, phreds , _ = read_fastq(fn, seqtype=String);
        noises = [mean(phred_to_p(i)) for i in phreds];
        #Calculate the expected proportion of sequences that are error free.
        expected_zero_errors = mean([pdf(Poisson(noises[i]*length(seqs[i])),0) for i in 1:length(seqs)]);
    end

    counts = countmap(seqs);
    seq_keys = keys(counts);
    seq_freqs_all = reverse(sort([(counts[k],k,kmer_count(k,6)) for k in seq_keys]));
    seq_freqs = seq_freqs_all[[k[1]>=2 for k in seq_freqs_all]];

    if length(seq_freqs) > 0
        current_set = [seq_freqs[1]];
        for i in 2:length(seq_freqs)
            dists = [corrected_kmer_dist_full(k[3],seq_freqs[i][3]) for k in current_set]
            neighbour_set = dists.<=neigh_thresh

            #if no neighbours within radius, add to current_set
            if sum(neighbour_set)==0
                push!(current_set,seq_freqs[i])
            else
                if method == 2
                    inds = collect(1:length(current_set))[neighbour_set]
                    best = indmax(current_set[inds])
                    p_val = (1-cdf(Poisson((current_set[inds][best][1]/expected_zero_errors)*err_rate),seq_freqs[i][1]))*length(seq_freqs[i][2])
                    if p_val < alpha
                        push!(current_set,seq_freqs[i])
                    end
                end
                #Using same as above. Will add in homopolymer test.
                if method == 3
                    #Need to look through all inds. For each, calculate p-value. Take the largest p-value. If that is <thresh, then keep cluster.
                    #This will use a regular error rate if HP==false, but a much larger error rate if HP==true.

                    inds = collect(1:length(current_set))[neighbour_set]
                    best = indmax(current_set[inds])
                    p_val = (1-cdf(Poisson((current_set[inds][best][1]/expected_zero_errors)*err_rate),seq_freqs[i][1]))*length(seq_freqs[i][2])
                    if p_val < alpha
                        push!(current_set,seq_freqs[i])
                    end
                end
            end
        end

        println("Using ",100*sum([k[1] for k in current_set])/length(seqs),"% of reads. If this is lower than 5%, use RAD instead.")

        frequencies = zeros(length(current_set))
        for s in seq_freqs_all
            dists = [corrected_kmer_dist_full(k[3],s[3],k=5) for k in current_set]
            frequencies[indmin(dists)] += s[1]
        end
        return [k[2] for k in current_set],frequencies
    else #catching the condition where there are no duplicates
        warn("Failing completely. Use RAD instead.")
        return [],[]
    end
end
