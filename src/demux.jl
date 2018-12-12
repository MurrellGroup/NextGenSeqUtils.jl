const IUPACnumDict = Dict{Char,Int8}(Dict())
IUPACnumDict['A']=1
IUPACnumDict['C']=2
IUPACnumDict['G']=3
IUPACnumDict['T']=4
IUPACnumDict['U']=5
IUPACnumDict['R']=6
IUPACnumDict['Y']=7
IUPACnumDict['S']=8
IUPACnumDict['W']=9
IUPACnumDict['K']=10
IUPACnumDict['M']=11
IUPACnumDict['B']=12
IUPACnumDict['D']=13
IUPACnumDict['H']=14
IUPACnumDict['V']=15
IUPACnumDict['N']=16

IUPACnucLookup=[[true false false false false true false false true false true false true true true true]
[false true false false false false true true false false true true false true true true]
[false false true false false true false true false true false true true false true true]
[false false false true true false true false true true false true true true false true]
[false false false true true false true false true true false true true true false true]
[true false true false false true false true true true true true true true true true]
[false true false true true false true true true true true true true true true true]
[false true true false false true true true false true true true true true true true]
[true false false true true true true false true true true true true true true true]
[false false true true true true true true true true false true true true true true]
[true true false false false true true true true false true true true true true true]
[false true true true true true true true true true true true true true true true]
[true false true true true true true true true true true true true true true true]
[true true false true true true true true true true true true true true true true]
[true true true false false true true true true true true true true true true true]
[true true true true true true true true true true true true true true true true]];

function toIUPACnum(strin::String)
    return [IUPACnumDict[i] for i in strin]
end

function IUPAC_equals(a1,a2)
    return IUPACnucLookup[a1,a2]
end

"""Needleman-Wunch, with IUPAC codes, and end gaps penalized slightly less."""
function IUPAC_nuc_edit_dist(s1::String, s2::String; edge_reduction = 0.9999)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 0.0

    
    s1numArr = toIUPACnum(s1)
    s2numArr = toIUPACnum(s2)
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1,:] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            elseif IUPAC_equals(s1numArr[i-1],s2numArr[j-1])
                diag = arr[i-1, j-1] + match_cost-0.00001
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = indmax(scores)
            if best == 1
                arr[i, j] = diag
            elseif best ==2
                arr[i, j] = del
            elseif best ==3
                arr[i, j] = ins
            end
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]
    return alignedScore
end

"""Needleman-Wunch, with IUPAC codes, and end gaps penalized slightly less."""
function IUPAC_nuc_nw_align(s1::String, s2::String; edge_reduction = 0.9999)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 0.0

    
    s1numArr = toIUPACnum(s1)
    s2numArr = toIUPACnum(s2)
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1,:] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            elseif IUPAC_equals(s1numArr[i-1],s2numArr[j-1])
                diag = arr[i-1, j-1] + match_cost-0.00001
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = indmax(scores)
            if best == 1
                arr[i, j] = diag
            elseif best ==2
                arr[i, j] = del
            elseif best ==3
                arr[i, j] = ins
            end
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]

    # return arr

    # the trace endpoing will need to be generalized if we want to
    # allow overhang.
    trI, trJ = length(s1arr), length(s2arr)

    # First compute the trace running backwards. Initialized to the
    # maximum size. With unit penalties, is it possible to tell how
    # big this should be in advance?
    backtrace = Array{Int}(length(s1arr) + length(s2arr))
    # If you hit any boundary, you have to run all the way to the side.
    btInd = 1
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = traceArr[trI, trJ]
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
        elseif backtrace[btInd] == 2
            trJ += -1
        elseif backtrace[btInd] == 3
            trI += -1
        end
        btInd += 1
    end

    # If you hit the boundaries not at the top left corner.
    while trI > 0
        backtrace[btInd] = 3
        btInd += 1
        trI += -1
    end

    while trJ > 0
        backtrace[btInd] = 2
        btInd += 1
        trJ += -1
    end

    backtrace = backtrace[1:btInd-1]

    # This will need to be generalized to work on non-strings. Not
    # important.
    ali1arr = Array{Char}(length(backtrace))
    ali2arr = Array{Char}(length(backtrace))
    ind1 = 1
    ind2 = 1
    for i in 1:length(backtrace)
        tr = backtrace[length(backtrace)-(i-1)]
        if tr == 1
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = s2arr[ind2]
            ind1 += 1
            ind2 += 1
        elseif tr == 2
            ali2arr[i] = s2arr[ind2]
            ali1arr[i] = '-'
            ind2 += 1
        elseif tr == 3
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = '-'
            ind1 += 1
        end
    end
    return join(ali1arr), join(ali2arr), alignedScore
end


function banded_edit_dist(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)
    # calculate band width
    avg_len = (length(s1) + length(s2)) / 2
    bandwidth = Int(ceil(band_coeff * sqrt(avg_len)))
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 1.0
    
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    
    # dimension difference added to band length on one side
    # to give padding for large sequence length differences
    dim_diff = length(s1arr) - length(s2arr)
    maxlen = max(length(s1), length(s2))
    arr = zeros(maxlen + 1, 2*bandwidth + 1 + abs(dim_diff))
    traceArr = zeros(Int, size(arr))
    for i in 1:bandwidth+1
        add_to_band!(arr, edge_reduction*del_cost*(i-1), i, 1, bandwidth, dim_diff)
        add_to_band!(arr, edge_reduction*ins_cost*(i-1), 1, i, bandwidth, dim_diff)
    end

    for i in 2:length(s1arr)+1
        # keep j within 2*bandwidth + 1 band of i, but keep bandwidth sized buffer on edges
        lo = 0
        # to keep a big enough padding around the edges when there's a big length mismatch
        if i <= bandwidth + 1 + max(0, dim_diff)
            lo = 2
        else
            lo = i - bandwidth - max(0, dim_diff)
        end
        for j in lo:length(s2arr)+1
            # keep j within band
            if !(in_band(i, j, bandwidth, dim_diff))
                break
            end
            
            if in_band(i-1, j-1, bandwidth, dim_diff) && s1arr[i-1] == s2arr[j-1]
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + match_cost
            else
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + mismatch_cost
            end
            
            # to handle the lower edge penalties.
            delMult = (i == length(s1arr)+1) ? edge_reduction : 1
            insMult = (j == length(s2arr)+1) ? edge_reduction : 1

            ins = in_band(i-1, j, bandwidth, dim_diff) ? 
                    (get_band_val(arr, i-1, j, bandwidth, dim_diff)+(ins_cost*insMult)) : NaN
            del = in_band(i, j-1, bandwidth, dim_diff) ? 
                    (get_band_val(arr, i, j-1, bandwidth, dim_diff)+(del_cost*delMult)) : NaN

            scores = [diag, del, ins]
            best = indmax(scores)
            add_to_band!(arr, scores[best], i, j, bandwidth, dim_diff)
            add_to_band!(traceArr, best, i-1, j-1, bandwidth, dim_diff)                                         
        end
    end
    alignedScore = get_band_val(arr, length(s1arr)+1, length(s2arr)+1, bandwidth, dim_diff)                                                                      
    
    return alignedScore
end


#Indmax in 2D, returning the coordinates
function maxcoord(mat)
    return ind2sub(mat,indmax(mat))
end

#THIS PRIMER MATCHING FUNCTION DOESN'T ASSUME WE KNOW PRIMER PAIRS
function matchPrimer(qseq::Tuple{Any,Any,Any},forwardPrimers,revPrimers,forwardNames,revNames)
    fwd_qseq = qseq;
    fwd_seq = fwd_qseq[1];
    rev_seq = reverse_complement(fwd_seq);
    rev_qseq = (rev_seq, fwd_qseq[2], fwd_qseq[3]); # The second and third element remain the same
    fwdSeqFwdPrim = [IUPAC_nuc_edit_dist(fwd_seq[1:length(prim)],prim) for prim in forwardPrimers];
    revSeqFwdPrim = [IUPAC_nuc_edit_dist(rev_seq[1:length(prim)],prim) for prim in forwardPrimers];
    fwdSeqRevPrim = [IUPAC_nuc_edit_dist(fwd_seq[1:length(prim)],prim) for prim in revPrimers];
    revSeqRevPrim = [IUPAC_nuc_edit_dist(rev_seq[1:length(prim)],prim) for prim in revPrimers];

    fwdDirectionScores = zeros(length(forwardPrimers),length(revPrimers))
    revDirectionScores = zeros(length(forwardPrimers),length(revPrimers))
    
    for i in 1:length(forwardPrimers)
        for j in 1:length(revPrimers)
            fwdDirectionScores[i,j]=fwdSeqFwdPrim[i]+revSeqRevPrim[j];
            revDirectionScores[i,j]=revSeqFwdPrim[i]+fwdSeqRevPrim[j];
        end
    end
    
    fwd,back = maximum(fwdDirectionScores),maximum(revDirectionScores)
    
    #fwdDirectionScores = fwdSeqFwdPrim.+revSeqRevPrim;
    #revDirectionScores = revSeqFwdPrim.+fwdSeqRevPrim;
    fwd,back = maximum(fwdDirectionScores),maximum(revDirectionScores)
    fwdDir = fwd>back
    if fwdDir
        primInd = maxcoord(fwdDirectionScores)
        worstScore = minimum([fwdSeqFwdPrim[primInd[1]],revSeqRevPrim[primInd[2]]])
        retseq = fwd_qseq
    else
        primInd = maxcoord(revDirectionScores)
        worstScore = minimum([revSeqFwdPrim[primInd[1]],fwdSeqRevPrim[primInd[2]]])
        retseq = rev_qseq
    end
    return primInd,forwardNames[primInd[1]]*"_"revNames[primInd[2]],retseq,round(worstScore),(!fwdDir)
end

function deMux(qseqs::Tuple{Array,Array,Array},forwardPrimers,reversePrimers, forwardNames,reverseNames)
    matched = [matchPrimer(qseq,forwardPrimers,reversePrimers,forwardNames,reverseNames) for qseq in zip(qseqs...)];
    return [[i[2],i[3], i[5], i[4]] for i in matched]
end

function demux_fastx(fwdPrimers, revPrimers; minLen = 2300, errorRate = 0.01,
    path = "/Users/vrkumar/Downloads/",
    filename = "LP97.fastq",
    mismatchThresh = 2
)
    forwardNames = [i[1] for i in fwdPrimers];
    reverseNames = [i[1] for i in revPrimers];
    forwardPrimers = [i[2] for i in fwdPrimers];
    reversePrimers = [i[2] for i in revPrimers];

    newName = path*filename*"_"*string(minLen)*"_"*string(errorRate)*".fastq"
    usearch_filter(path*filename,newName, errorRate=errorRate,minLength=minLen,labelPrefix="seq",errorOut = true)

    if filename[end] == 'a'
        names, seqs = read_fasta_with_names(newName)
        phreds = nothing;
    else
        seqs, phreds, names = read_fastq(newName);
    end

    deMuxed = deMux(seqs,forwardPrimers,reversePrimers,forwardNames,reverseNames);

    primerPairs = [i[1] for i in deMuxed];
    primerMatchedSeqs = [i[2] for i in deMuxed];
    reversed = [i[3] for i in deMuxed];
    scores = [i[4] for i in deMuxed]
    
    if phreds != nothing
        out_phreds = [reversed[i] ? reverse(phreds[i]) : phreds[i] for i in 1:length(reversed)]
    else
        out_phreds = [0 for i in 1:length(reversed)]
    end
        
    out = Dict()
    for i in 1:length(primerPairs)
        if scores[i] > mismatchThresh
            continue
        end
        if !(haskey(out, primerPairs[i]))
            out[primerPairs[i]] = [[],[],[]]
        end
        push!(out[primerPairs[i]][1], primerMatchedSeqs[i])
        push!(out[primerPairs[i]][2], out_phreds[i])
        push!(out[primerPairs[i]][3], names[i])
    end
    
    for p in keys(out)
        if phreds != nothing
            write_fastq(path*p*".fastq", out[p][1], Array{Array{Int8,1},1}(out[p][2]), names = out[p][3])
        else
            write_fasta(path*p*".fastq", out[p][1], names = out[p][3])
        end
    end
    
end
