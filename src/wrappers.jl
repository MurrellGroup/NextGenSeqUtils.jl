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

function substring_coords(query,ref)
    start_coord = searchindex(ref,query)
    if start_coord == 0
        return false        
    end
    return (start_coord,start_coord+length(query)-1)
end

function trim_fastq_entry(trimmed_seq,fastq_read,phred_arr)
    trimmed_seq = uppercase(trimmed_seq)
    fastq_read = uppercase(fastq_read)
    #first look for seq in current orientation
    revcomp = false
    try_coords = substring_coords(trimmed_seq,fastq_read)
    #then try the reverse complement
    if try_coords==false
        revcomp = true
        fastq_read = reverse_complement(fastq_read)
        try_coords = substring_coords(trimmed_seq,fastq_read)
        if try_coords==false
            return false
        else
            phred_arr = reverse(phred_arr)
        end
    end
    return fastq_read[try_coords[1]:try_coords[2]],phred_arr[try_coords[1]:try_coords[2]]
end

function array2index_dict(arr)
    dic = Dict{String,Int64}()
    for i in 1:length(arr)
        dic[arr[i]]=i
    end
    return dic
end

function usearch_trim_fastq_with_phreds(inFastqPath, referencePath; divergence=0.9, outFastqPath="")
    #READ IN THE FASTQ FILE
    fastqseqs, fastqphreds, fastqnames = read_fastq(inFastqPath)
    fastq_name_dic = array2index_dict(fastqnames)
    
    #USEARCH TO TRIM, BUT ONLY STORES .FASTA OUTPUT
    tempOutPath = "/tmp/usearch9tempOutFile"
    run(`$(PATHS.usearch) -usearch_global $inFastqPath -db $referencePath -id $divergence -fastapairs $tempOutPath -strand both -top_hit_only`)
    names, seqs = read_fasta_with_names(tempOutPath)
    ref_names = [names[i*2] for i in 1:Int64(round(length(names)/2))]
    
    names2nd = [names[i*2-1] for i in 1:Int64(round(length(names)/2))]
    seqs2nd = [degap(seqs[i*2-1]) for i in 1:Int64(round(length(names)/2))]
    
    #trim_fastq_entry returns read,phred
    trimmed_fastq = [(fastqnames[fastq_name_dic[names2nd[i]]],trim_fastq_entry(seqs2nd[i],fastqseqs[fastq_name_dic[names2nd[i]]],fastqphreds[fastq_name_dic[names2nd[i]]])) for i in 1:length(names2nd)]

    if outFastqPath != ""
        write_fastq(outFastqPath, [i[2][1] for i in trimmed_fastq], [i[2][2] for i in trimmed_fastq],names=[i[1] for i in trimmed_fastq])
    end
    
    return ref_names,trimmed_fastq
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
function mafft_consensus(seqs::Vector{T},; kwargs...) where T <: BioSequence
    mktempdir() do mydir
        seqfile = string(mydir, "/sequences.fasta")
        mafftout = string(mydir, "/mafft.fasta")

        write_fasta(string(mydir, "/sequences.fasta"), seqs)
        mafft(seqfile, mafftout; kwargs...)
        aligned = Array{String}(read_fasta(mafftout))
        return alignment_consensus(aligned)
    end
end
