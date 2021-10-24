# copypasta from kemal for fasta and fastq i/o

# i/o functions use file parsing from BioSequences package

### FASTA files

"""
    read_fasta_records(filename)

Read .fasta file contents.
"""
function read_fasta_records(filename)
    stream = open(FASTA.Reader, filename)
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    close(stream)
    return records
end

"""
    read_fasta(filename; seqtype=String)

Read .fasta file contents, parse, and return names and sequences as type `seqtype`.
"""
function read_fasta(filename; seqtype=String)
    records = read_fasta_records(filename)
    return [FASTA.identifier(r) for r in records], seqtype[FASTA.sequence(seqtype, r) for r in records]
end

#for legacy code
const read_fasta_with_names = read_fasta

"""
    read_fasta_with_names_in_other_order(filename; seqtype=String)

Read .fasta file contents, parse, and return sequences as type `seqtype` with names.
Now that's what I call convenience.
"""
function read_fasta_with_names_in_other_order(filename::String; seqtype=String)
    names, seqs = read_fasta_with_names(filename, seqtype=seqtype)
    return seqs, names
end

"""
    write_fasta(filename::String, seqs; names = String[])

Write given `seqs` and optional `names` to a .fasta file with given filepath.
"""
function write_fasta(filename::String, seqs; names = String[])
    if length(names) > 0 && length(names) != length(seqs)
        error("number of sequences does not match number of names")
    end
    if length(names) == 0
        names = ["seq_$i" for i in 1:length(seqs)]
    end
    stream = open(FASTA.Writer, filename)
    for (name, seq) in zip(names, seqs)
        write(stream, FASTA.Record(name, seq))
    end
    close(stream)
end


### FASTQ files

"""
    read_fastq_records(filename)

Read .fastq file contents.
"""
function read_fastq_records(filename)
    stream = open(FASTQ.Reader, filename)
    records = FASTQ.Record[]
    for record in stream
        if any([q < 0 for q in record.quality])
            error("$(record.name) in $filename contains negative phred values")
        end
        push!(records, record)
    end
    close(stream)
    return records
end

"""
    read_fastq(filename; seqtype=String)

Read .fastq file contents, parse, and return sequences as `seqtype` type, phreds, and names.
"""
function read_fastq(filename; seqtype=String, min_length=nothing, max_length=nothing, err_rate=nothing)
    records = read_fastq_records(filename)
    seqs = seqtype[]
    phreds = Vector{Phred}[]
    names = String[]
    for record in records
		if err_rate != nothing && mean(phred_to_p(Array{Int8,1}(FASTQ.quality(record, :sanger)))) > err_rate
			continue
		end
		if min_length != nothing && length(FASTQ.sequence(seqtype, record)) < min_length
			continue
		end
		if max_length != nothing && length(FASTQ.sequence(seqtype, record)) > max_length
			continue
        end
		push!(seqs, FASTQ.sequence(seqtype, record))
        push!(phreds, FASTQ.quality(record, :sanger))
        push!(names, FASTQ.identifier(record))
    end
    return seqs, phreds, names
end

"""
    write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                     names=String[], LongSequence = false)

Write given sequences, phreds, names to .fastq file with given file path.
If `names` not provided, gives names 'seq_1', etc.
"""
function write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                     names=String[], LongSequence = false)
    if !LongSequence
        seqs = [LongCharSeq(s) for s in seqs]
    end
    stream = open(FASTQ.Writer, filename)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, q, n) in zip(seqs, phreds, names)
        i += 1
        write(stream, FASTQ.Record(n, s, q))
    end
    close(stream)
end

#------Chunked IO funcions--------

"""
function chunked_fastq_apply(fpath, func::Function; chunk_size=10000, f_kwargs = [], verbose = false)

Pass a function of the form `func(seqs::Array{Any,1}, phreds::Array{Phred,1}, names::Array{Any,1})` to apply to a large FASTQ file chunk by chunk.
The input file can be Gzipped or uncompressed. Additional kwargs can be passed to the function via f_kwargs.
"""
function chunked_fastq_apply(fpath, func::Function; chunk_size=10000, f_kwargs = [], verbose = false)
	if endswith(fpath, ".gz")
    	reader = FASTQ.Reader(GzipDecompressorStream(open(fpath)))
	else
		reader = FASTQ.Reader(open(fpath))
	end
	if !hasmethod(func, Tuple{Array{Any,1}, Array{Phred,1}, Array{Any,1}})
		@error "Function argument must accept func(seqs::Array{Any,1}, phreds::Array{Phred,1}, names::Array{Any,1})!"
	end
    seqs, phreds, names = [], Vector{Phred}[], []
    i = 0
    record = FASTQ.Record()
	results = []
    while !eof(reader)
        read!(reader, record)
        push!(seqs, FASTQ.sequence(String, record))
        push!(phreds, FASTQ.quality(record, :sanger))
        push!(names, FASTQ.identifier(record))
        i += 1
        if i == chunk_size
            #apply func...
            res = func(seqs, phreds, names; f_kwargs...)
			push!(results, res)
            seqs, phreds, names = [], Vector{Phred}[], []
            i = 0
            if verbose print(".") end
        end
    end
    close(reader)
    if i > 0
        #apply func...
        res = func(seqs, phreds, names; f_kwargs...)
		push!(results, res)
    end
    if verbose print(".") end
	return results
end
