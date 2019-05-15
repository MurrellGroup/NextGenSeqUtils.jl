#-------PHREDS AND PROBS--------
#copypasta from rifraf

const Phred = Int8
const Prob = Float64
const LogProb = Float64

const MIN_PHRED = Phred(1)
const MAX_PHRED = Phred(Int('~') - 33)

"""
    phred_to_log_p(x)

Conversion from phred value to log probability value.
"""
@generated function phred_to_log_p(x)
    return quote
        return x / (-10.0)
    end
end

"""
    phred_to_p(q::Phred)

Conversion from phred value to real probability value.
"""
function phred_to_p(q::Phred)
    return exp10(phred_to_log_p(q))
end

"""
    phred_to_p(x::Vector{Phred})

Conversion from phred values to real probability values in vector.
"""
function phred_to_p(x::Vector{Phred})
    return exp10.(phred_to_log_p(x))
end

"""
    p_to_phred(p::Prob)

Conversion from real probability value to phred value.
"""
function p_to_phred(p::Prob)
    return Phred(min(round(-10.0 * log10(p)), MAX_PHRED))
end

"""
    p_to_phred(p::Prob)

Conversion from real probability values to phred values in vector.
"""
function p_to_phred(x::Vector{LogProb})
    return Phred[p_to_phred(p) for p in x]
end

"""
    error_probs_to_phreds(eps::Vector{Float64})

Conversion from error probability values to phred values in vector.
"""
function error_probs_to_phreds(eps::Vector{Float64})
    return [Phred(min(round(-10.0 * log10(p)),99)) for p in eps]
end

#-------ALTERNATE FILTER FUNCTION----------

"""
    quality_filter(infile, outfile=join(split(infile, ".")[1:end-1], ".") * ".filt.fasta"; 
                   errorRate = 0.01, minLength = 100, labelPrefix = "seq", errorOut = true)

Writes file with sequences from input file that have all sites within error rate margin.
"""
function quality_filter(infile, outfile=join(split(infile, ".")[1:end-1], ".") * ".filt.fasta" ; 
                        errorRate = 0.01, minLength = 100, labelPrefix = "seq", errorOut = true)
    seqs, scores, names = read_fastq(infile)

	p_vals = [phred_to_p(score) for score in scores]
    
    inds = quality_filter_inds(seqs, scores, errorRate=errorRate, minLength=minLength)
        
    #names = names[inds]
    if errorOut == true
        names = ["$labelPrefix$(i)|ee=$(mean(p_vals[i]))" for i in inds]
    else
        names = ["$labelPrefix$(i)" for i in 1:length(inds)]
    end
    
    if last(outfile) == 'a'
        write_fasta(outfile, seqs[inds], names=names)
    else
        write_fastq(outfile, seqs[inds], scores[inds], names=names)
    end
end

"""
    quality_filter(seqs::Array{String, 1}, scores::Array, names::Array{String, 1}; errorRate=0.01, minLength=0)

Returns sub arrays with sequences that have site-wise error all within given margin.
"""
function quality_filter(seqs::Array{String, 1}, scores::Array, names::Array{String, 1}; errorRate=0.01,minLength=0)
    inds = quality_filter_inds(seqs, scores, errorRate=errorRate, minLength=minLength)
    return seqs[inds], scores[inds], names[inds]
end

"""
    quality_filter_inds(seqs::Array{String, 1}, scores::Array; errorRate=0.01,minLength=100)

Returns indices of sequences with site-wise error within given margin.
"""
function quality_filter_inds(seqs::Array{String, 1}, scores::Array; errorRate=0.01,minLength=100)
    p_vals = [phred_to_p(score) for score in scores]
    inds = filter!(x -> (mean(p_vals[x]) < errorRate && (length(seqs[x]) >= minLength)) , collect(1:length(seqs)))
    return inds
end

# temporary: for compatibility
seq_filter = quality_filter

#--------QUALITY INTERROGATING FUNCTIONS, GIVEN EITHER EE NAME TAGGING OR FASTQ FORMAT ------- 

"""
    length_vs_qual(fasta_path; plot_title = "Length Vs Errors")

Creates plot of lengths of sequences vs. mean error rates of sequences.
"""
function length_vs_qual(fasta_path; plot_title = "Length Vs Errors", alpha=0.3)
	
	if last(fasta_path) == 'a'
		names, seqs = read_fasta_with_names(fasta_path)
    	error_rates = [parse(Float64, split(i,['|', '='])[3]) for i in names]
	else
		seqs, errs, names = read_fastq(fasta_path)
		error_rates = [mean(phred_to_p(e)) for e in errs]
    end

	# To remove that annoying rounding that happens in usearch
    # dithered_rates = [(rand()*10-5)*((i%1)==0.0&&i>90.0)+i for i in error_rates]
	# Following line commented since usearch is no longer used    
	# dithered_rates = [(1+rand()*0.1-0.05)*i for i in error_rates]
    string_lengths = [length(i) for i in seqs]
    plot(string_lengths, error_rates, ".", alpha=alpha)
    title(plot_title)
    xlabel("Length")
    ylabel("Predicted Error Rate")
    # return dithered_rates
end

"""
    qual_hist(fasta_path; plot_title = "Quality Hist")

Creates histogram of mean quality scores per sequence.
"""
function qual_hist(fasta_path; plot_title = "Quality Hist")
	if last(fasta_path) == 'a'
		names, seqs = read_fasta_with_names(fasta_path)
    	error_rates = [parse(Float64, split(i,['|', '='])[3]) for i in names]
	else
		seqs, errs, names = read_fastq(fasta_path)
		error_rates = [mean(phred_to_p(e)) for e in errs]
    end
    # to remove that annoying rounding that happens in usearch
    # Following line commented since usearch is no longer used
	#dithered_rates = [(rand()*10-5)*((i%1)==0.0&&i>90.0)+i for i in error_rates]
    string_lengths = [length(i) for i in seqs]
    plt[:hist](error_rates, 40)
    title(plot_title)
    xlabel("Error Rate")
    ylabel("Frequency")
    # return dithered_rates
end

export fastq_filter
function fastq_filter(infile, outfile; 
        error_rate = 0.01, min_length = 30, max_length = 1000000,
        label_prefix = "seq", error_out = true)
    
    println("Reading file...")
    seqs, scores, names = read_fastq(infile)
    
    println("Calculating filter...")
    lengths = length.(seqs)
    mean_errors = [mean(NextGenSeqUtils.phred_to_p.(score)) for score in scores]
    inds = [1:length(seqs);][(lengths .< max_length) .& (lengths .> min_length) .& (mean_errors .< error_rate)]
    
    if error_out == true
        names = ["$label_prefix$(i)|ee=$(mean_errors[i])" for i in inds]
    else
        names = ["$label_prefix$(i)" for i in 1:length(inds)]
    end
    
    println("Writing file...")
    if last(outfile) == 'a'
        write_fasta(outfile, seqs[inds], names=names)
    else
        write_fastq(outfile, seqs[inds], scores[inds], names=names)
    end
end
