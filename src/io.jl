# copypasta from kemal for fasta and fastq i/o

### FASTA files

function read_fasta_records(filename)
    stream = open(FASTA.Reader, filename)
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

function read_fasta(filename; seqtype=String)
    records = read_fasta_records(filename)
    return seqtype[FASTA.sequence(seqtype, r) for r in records]
end

function read_fasta_with_names(filename; seqtype=String)
    records = read_fasta_records(filename)
    return [FASTA.identifier(r) for r in records], seqtype[FASTA.sequence(seqtype, r) for r in records]
end

function read_fasta_with_names_in_other_order(filename::String; seqtype=String)
    names, seqs = read_fasta_with_names(filename, seqtype=seqtype)
    return seqs, names
end

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

function read_fastq_records(filename)
    stream = open(FASTQ.Reader, filename)
    records = FASTQ.Record[]
    for record in stream
        if any([q < 0 for q in record.quality])
            error("$(record.name) in $filename contains negative phred values")
        end
        push!(records, record)
    end
    return records
end

function read_fastq(filename; seqtype=String)
    records = read_fastq_records(filename)
    seqs = seqtype[]
    phreds = Vector{Phred}[]
    names = String[]
    for record in records
        push!(seqs, FASTQ.sequence(seqtype, record))
        push!(phreds, FASTQ.quality(record, :sanger))
        push!(names, FASTQ.identifier(record))
    end
    return seqs, phreds, names
end

function write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                     names=String[], DNASeqType = false)
    if !DNASeqType
        seqs = [DNASequence(s) for s in seqs]
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
