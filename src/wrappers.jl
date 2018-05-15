#------USEARCH WRAPPERS---------

"""
    usearch_filter(inFastqPath, outFastXPath;
                   errorRate=0.01, minLength=100,
                   labelPrefix="seq", errorOut=true)

Julia wrapper for usearch error rate filter. Input file may be .fastq or .fasta.
"""
function usearch_filter(inFastqPath, outFastXPath;
                        errorRate=0.01, minLength=100,
                        labelPrefix="seq", errorOut=true)
    if last(outFastXPath) == 'q'
        outputTypeString = "-fastqout"
    else
        if last(outFastXPath)=='a'
            outputTypeString = "-fastaout"
            else error("Output filepath must end in 'q' (.fastq) or 'a' (.fasta)")
        end
    end
    errorOutString = ""
    if errorOut == true
        errorOutString = "-fastq_eeout"
    end
    run(`$(PATHS.usearch) -fastq_filter $inFastqPath -fastq_maxee_rate $errorRate -fastq_qmax 99 -fastq_minlen $minLength $outputTypeString $outFastXPath -relabel $labelPrefix $errorOutString`)
end


"""
    usearch_ref(inFastXPath, outFastaPath, referenceString;
                     divergence=0.9)

Julia wrapper for usearch.
"""
function usearch_ref(inFastXPath, outFastaPath, referenceString;
                     divergence=0.9)
    tempRefPath = "/tmp/usearch9tempRefFile"
    tempOutPath = "/tmp/usearch9tempOutFile"
    write_fasta(tempRefPath,[DNASequence(referenceString)])
    run(`$(PATHS.usearch) -usearch_global $inFastXPath -db $tempRefPath -id $divergence -fastapairs $tempOutPath -strand both`)
    names, seqs = read_fasta_with_names(tempOutPath)
    names2nd = [names[i * 2 - 1] for i in 1:Int64(round(length(names) / 2))]
    seqs2nd = [degap(seqs[i * 2 - 1]) for i in 1:Int64(round(length(names) / 2))]
    write_fasta(outFastaPath, seqs2nd, names=names2nd)
    return names2nd, seqs2nd
end

"""
    usearch_ref_database(inFastXPath, outFastaPath, referencePath;
                         divergence=0.9)

Julia wrapper for usearch.
"""
function usearch_ref_database(inFastXPath, outFastaPath, referencePath;
                              divergence=0.9)
    tempOutPath = "/tmp/usearch9tempOutFile"
    run(`$(PATHS.usearch) -usearch_global $inFastXPath -db $referencePath -id $divergence -fastapairs $tempOutPath -strand both -top_hit_only`)
    names, seqs = read_fasta_with_names(tempOutPath)
    names2nd = [names[i*2-1] for i in 1:Int64(round(length(names)/2))]
    seqs2nd = [degap(seqs[i*2-1]) for i in 1:Int64(round(length(names)/2))]
    write_fasta(outFastaPath, seqs2nd, names=names2nd)
    return names2nd, seqs2nd
end

#-------MAFFT WRAPPER--------
"""
    mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)

Julia wrapper for mafft.
"""
function mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
    mafft = length(path) == 0 ? PATHS.mafft : path
    flagstrings = []
    for flag in flags
        push!(flagstrings, "--$flag")
    end
    args = []
    for (k, v) in kwargs
        push!(args, "--$k")
        push!(args, "$v")
    end
    run(`$mafft $flagstrings $args --quiet --adjustdirection --progress /tmp/mafft.progress --out $outpath $inpath`)
end

"""
    mafft_consensus{T<:BioSequence}(seqs::Vector{T}; kwargs...)

Julia wrapper for mafft.
"""
function mafft_consensus{T<:BioSequence}(seqs::Vector{T}; kwargs...)
    mktempdir() do mydir
        seqfile = string(mydir, "/sequences.fasta")
        mafftout = string(mydir, "/mafft.fasta")

        write_fasta(string(mydir, "/sequences.fasta"), seqs)
        mafft(seqfile, mafftout; kwargs...)
        aligned = Array{String}(read_fasta(mafftout))
        return alignment_consensus(aligned)
    end
end
