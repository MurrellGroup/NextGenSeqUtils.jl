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
            best = findmax(scores)[2]
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
            best = findmax(scores)[2]
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
    backtrace = Array{Int}(undef, length(s1arr) + length(s2arr))
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
    ali1arr = Array{Char}(undef, length(backtrace))
    ali2arr = Array{Char}(undef, length(backtrace))
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
            best = findmax(scores)[2]
            add_to_band!(arr, scores[best], i, j, bandwidth, dim_diff)
            add_to_band!(traceArr, best, i-1, j-1, bandwidth, dim_diff)                                         
        end
    end
    alignedScore = get_band_val(arr, length(s1arr)+1, length(s2arr)+1, bandwidth, dim_diff)                                                                      
    
    return alignedScore
end


#Indmax in 2D, returning the coordinates
function maxcoord(mat)
    return findmax(mat)[2]
end

#THIS PRIMER MATCHING FUNCTION DOESN'T ASSUME WE KNOW PRIMER PAIRS
function matchPrimer(seq,forwardPrimers,revPrimers,forwardNames,revNames; epsilon = 2, trim = false, 
        forwardBarcodeMasks = nothing, reverseBarcodeMasks = nothing, barcoded = false)
    fwdseq = seq;
    revseq = reverse_complement(fwdseq);
    fwdSeqFwdPrim = [IUPAC_nuc_edit_dist(fwdseq[1:(length(prim) + epsilon)],prim, edge_reduction = 0.5) for prim in forwardPrimers];
    revSeqFwdPrim = [IUPAC_nuc_edit_dist(revseq[1:(length(prim) + epsilon)],prim, edge_reduction = 0.5) for prim in forwardPrimers];
    fwdSeqRevPrim = [IUPAC_nuc_edit_dist(fwdseq[1:(length(prim) + epsilon)],prim, edge_reduction = 0.5) for prim in revPrimers];
    revSeqRevPrim = [IUPAC_nuc_edit_dist(revseq[1:(length(prim) + epsilon)],prim, edge_reduction = 0.5) for prim in revPrimers];

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
    
    trunc_inds = nothing
    fwalign, revalign = nothing, nothing
    primInd = nothing
    
    if fwdDir
        primInd = maxcoord(fwdDirectionScores)
        worstScore = minimum([fwdSeqFwdPrim[primInd[1]],revSeqRevPrim[primInd[2]]])
        if trim 
            prim1 = forwardPrimers[primInd[1]]
            prim2 = revPrimers[primInd[2]]
            fwalign = IUPAC_nuc_nw_align(fwdseq[1:(length(prim1) + epsilon)],prim1, edge_reduction = 0.5)
            revalign = IUPAC_nuc_nw_align(revseq[1:(length(prim2) + epsilon)],prim2, edge_reduction = 0.5)
            ind1 = (findlast(x -> x != '-', fwalign[2])) + 1
            ind2 = length(fwdseq) - (findlast(x -> x != '-', revalign[2]))
            retseq = fwdseq[ind1:ind2]
            trunc_inds = ind1:ind2
        else
            retseq = fwdseq
        end
    else
        primInd = maxcoord(revDirectionScores)
        worstScore = minimum([revSeqFwdPrim[primInd[1]],fwdSeqRevPrim[primInd[2]]])
        if trim 
            prim1 = revPrimers[primInd[1]]
            prim2 = forwardPrimers[primInd[2]]
            fwalign = IUPAC_nuc_nw_align(revseq[1:(length(prim1) + epsilon)],prim1, edge_reduction = 0.5)
            revalign = IUPAC_nuc_nw_align(fwdseq[1:(length(prim2) + epsilon)],prim2, edge_reduction = 0.5)
            ind1 = (findlast(x -> x != '-', fwalign[2])) + 1
            ind2 = length(revseq) - (findlast(x -> x != '-', revalign[2]))
            retseq = revseq[ind1:ind2]
            trunc_inds = ind1:ind2
        else
            retseq = revseq
        end
    end
    
    barcode_score = 0
    if barcoded
        ptr = 1
        index = 0
        last_ptr = (findlast(x -> x != '-', fwalign[2]))
        while ptr <= last_ptr
            if fwalign[2][ptr] != '-' 
                index += 1
            end
            if forwardBarcodeMasks[primInd[1]][max(index,1)] && fwalign[1] != fwalign[2]
                barcode_score += 1
            end
            ptr += 1
        end
        
        ptr = 1
        index = 0
        last_ptr = (findlast(x -> x != '-', revalign[2]))
        while ptr <= last_ptr
            if reverseBarcodeMasks[primInd[2]][max(index,1)] && revalign[1] != revalign[2]
                barcode_score += 1
            end
            if revalign[2][ptr] != '-' 
                index += 1
            end
            ptr += 1
        end
    end
    
    return primInd,forwardNames[primInd[1]]*"_"*revNames[primInd[2]],retseq,round(worstScore),(!fwdDir), trunc_inds, barcode_score 
end

function deMux(seqs,forwardPrimers,reversePrimers, forwardNames,reverseNames; epsilon = 0, trim = false, return_all = false, forwardBarcodeMasks = nothing, reverseBarcodeMasks = nothing, barcoded=false)
    matched = [matchPrimer(seq,forwardPrimers,reversePrimers,forwardNames,reverseNames, epsilon = epsilon, trim = trim, forwardBarcodeMasks = forwardBarcodeMasks, reverseBarcodeMasks = reverseBarcodeMasks, barcoded=barcoded) for seq in seqs];
    if return_all
        return matched
    else
        return [[i[2],i[3], i[5], i[4]] for i in matched]
    end
end

function demux_fastx(fwdPrimers, revPrimers; minLen = 2300, errorRate = 0.01,
    path = "/Users/vrkumar/Downloads/",
    filename = "LP97.fastq",
    mismatchThresh = 4, epsilon = 2, barcodeThresh = 1,
    trim = true, forwardBarcodeMasks = nothing, reverseBarcodeMasks = nothing, barcoded=false 
)
    forwardNames = [i[1] for i in fwdPrimers];
    reverseNames = [i[1] for i in revPrimers];
    forwardPrimers = [i[2] for i in fwdPrimers];
    reversePrimers = [i[2] for i in revPrimers];

    newName = path*filename
    #usearch_filter(path*filename,newName, errorRate=errorRate,minLength=minLen,labelPrefix="seq",errorOut = true)
    
    if filename[end] == 'a'
        names, seqs = read_fasta_with_names(newName)
        #phreds = nothing;
    else
        seqs, phreds, names = read_fastq(newName);
    end
    
    if phreds != nothing
        seqs, phreds, names = quality_filter(seqs, phreds, names; errorRate=errorRate)
    end
    #seqs, phreds, names = length_filter(seqs, phreds, names, minLength=minLen)
    
    deMuxed = [matchPrimer(s,forwardPrimers,reversePrimers,forwardNames,reverseNames, epsilon=epsilon, trim=trim, forwardBarcodeMasks = forwardBarcodeMasks, reverseBarcodeMasks = reverseBarcodeMasks, barcoded=barcoded) for s in seqs];
    
    primers = [i[1] for i in deMuxed];
    primerPairs = [i[2] for i in deMuxed];
        
    primerMatchedSeqs = [i[3] for i in deMuxed];
    reversed = [i[5] for i in deMuxed];
    scores = [i[4] for i in deMuxed]
    trunc_inds = [i[6] for i in deMuxed]
    barcode_scores = [i[7] for i in deMuxed]
    
    if phreds != nothing
        out_phreds = [reversed[i] ? reverse(phreds[i]) : phreds[i] for i in 1:length(reversed)]
    else
        out_phreds = [0 for i in 1:length(reversed)]
    end
        
    out = Dict()
    for i in 1:length(primerPairs)
        println(scores[i], barcode_scores[i])
        if scores[i] < -mismatchThresh || barcode_scores[i] < -barcodeThresh
            continue
        end
        if !(haskey(out, primerPairs[i]))
            out[primerPairs[i]] = [[],[],[]]
        end
        push!(out[primerPairs[i]][1], primerMatchedSeqs[i])
        push!(out[primerPairs[i]][2], out_phreds[i][trunc_inds[i]])
        push!(out[primerPairs[i]][3], names[i])
    end
     
    for p in keys(out)
        if phreds != nothing
            println(length(out[p][1]))
            println(length(out[p][2]))
            write_fastq(path*p*".fastq", out[p][1], Array{Array{Int8,1},1}(out[p][2]), names = out[p][3])
        else
            write_fasta(path*p*".fastq", out[p][1], names = out[p][3])
        end
    end
    
end




#Dictionary-based matching
function ambig_positions(seq::String)
    return [1:length(seq);][[!(s in ['A','C','G','T']) for s in seq]]
end
function ambig_positions(seq::Array{Char})
    return ambig_positions(join(seq))
end
function IUPAC_char_equals(c1,c2)
    return IUPAC_equals(IUPACnumDict[c1],IUPACnumDict[c2])
end
function ambig_options(c::Char)
    return ['A','C','G','T'][[IUPAC_char_equals('A',c),IUPAC_char_equals('C',c),IUPAC_char_equals('G',c),IUPAC_char_equals('T',c)]]
end
function replace_positions(seq,pos,chars)
    chararr = collect(seq)
    for i in 1:length(pos)
        chararr[pos[i]] = chars[i]
    end
    return join(chararr)
end
function ambig_expand(seq)
    inds = ambig_positions(seq)
    options = [ambig_options(seq[i]) for i in inds]
    collec = [[]]
    for to_add in options
        collec = vcat([[vcat(el,[addit]) for addit in to_add] for el in collec]...)
    end
    return [replace_positions(seq,inds,coll) for coll in collec]
end
function all_mutants_of_all_ambigs(seq)
    ambig_seqs = ambig_expand(seq)
    return vcat([all_mutants(s) for s in ambig_seqs]...)
end
function mutate_arr_at_pos!(arr,pos,elem)
    arr[pos] = elem
    return arr
end
function all_mutants(seq)
    arr_seq = collect(seq)
    delsA = join.([vcat(deleteat!(copy(arr_seq),i),['A']) for i in 1:length(arr_seq)])
    delsC = join.([vcat(deleteat!(copy(arr_seq),i),['C']) for i in 1:length(arr_seq)])
    delsG = join.([vcat(deleteat!(copy(arr_seq),i),['G']) for i in 1:length(arr_seq)])
    delsT = join.([vcat(deleteat!(copy(arr_seq),i),['T']) for i in 1:length(arr_seq)])
    insertsA = join.([insert!(copy(arr_seq),i,'A')[1:end-1] for i in 1:length(arr_seq)])
    insertsC = join.([insert!(copy(arr_seq),i,'C')[1:end-1] for i in 1:length(arr_seq)])
    insertsG = join.([insert!(copy(arr_seq),i,'G')[1:end-1] for i in 1:length(arr_seq)])
    insertsT = join.([insert!(copy(arr_seq),i,'T')[1:end-1] for i in 1:length(arr_seq)])
    mutatesA = join.([mutate_arr_at_pos!(copy(arr_seq),i,'A')[1:end] for i in 1:length(arr_seq)])
    mutatesC = join.([mutate_arr_at_pos!(copy(arr_seq),i,'C')[1:end] for i in 1:length(arr_seq)])
    mutatesG = join.([mutate_arr_at_pos!(copy(arr_seq),i,'G')[1:end] for i in 1:length(arr_seq)])
    mutatesT = join.([mutate_arr_at_pos!(copy(arr_seq),i,'T')[1:end] for i in 1:length(arr_seq)])
    return vcat([seq],delsA,delsC,delsG,delsT,insertsA,insertsC,insertsG,insertsT,mutatesA,mutatesC,mutatesG,mutatesT)
end

export fast_primer_match
function fast_primer_match(seqs,primers; tol_one_error = true)
    #Only tolerates a single bp difference between primer and seq
    #Assumes all filtering primers are the same length.
    #Note this doesn't mean that the primers actually had to be the same length! 
    l = length(primers[1])
    matches = zeros(Int,length(seqs))
    noisy_primer_map = Dict{String,Int}()
    for (i,pr) in enumerate(primers)
        if tol_one_error
            targets = all_mutants_of_all_ambigs(pr)
        else
            targets = ambig_expand(pr)
        end
        for npr in targets
            noisy_primer_map[npr] = i
        end
    end
    for (i,s) in enumerate(seqs)
        if matches[i] == 0
            matches[i] = get(noisy_primer_map,s[1:l],0)
        end
        if matches[i] == 0
            matches[i] = -get(noisy_primer_map,reverse_complement(s[end-l+1:end]),0)
        end
    end
    return matches #Return sign is the rev_comp direction. Magnitude is the primer index.
end

export fast_primer_pair_match
function fast_primer_pair_match(seqs,fwd_primers,rev_primers; tol_one_error = true)
    fwd_matches = fast_primer_match(seqs,fwd_primers,tol_one_error=tol_one_error)
    rev_matches = fast_primer_match(seqs,rev_primers,tol_one_error=tol_one_error)
    problem_count = sum(fwd_matches.*rev_matches .> 0)
    if problem_count>0
        @warn "Inconsistencies: $(problem_count)"
    end
    keepers = fwd_matches.*rev_matches .< 0
    rev_comp_bool = fwd_matches .< 0
    return keepers,abs.(fwd_matches),abs.(rev_matches),rev_comp_bool
end

export demux_dict
function demux_dict(seqs,fwd_primers,rev_primers; verbose = true, phreds = nothing, tol_one_error = true)
    if rev_primers == nothing
        fwd_matches = fast_primer_match(seqs,fwd_primers,tol_one_error=tol_one_error)
        rev_comp_bool = fwd_matches .< 0
        keepers = abs.(fwd_matches) .> 0
        fwd_matches = abs.(fwd_matches)
        pair_keeps = fwd_matches[keepers]        
    else
         keepers,fwd_matches,rev_matches,rev_comp_bool = fast_primer_pair_match(seqs,fwd_primers,rev_primers,tol_one_error=tol_one_error)
        f_keeps = fwd_matches[keepers]
        r_keeps = rev_matches[keepers]
        pair_keeps = [(f_keeps[i],r_keeps[i]) for i in 1:length(f_keeps)]
    end
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end
    
    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            end
        end
        return seq_dict
    end
end

export primer_trim
function primer_trim(seq,primer; buffer = 3)
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq[1:length(primer)+buffer])
    gapbool = reverse(collect(a1)) .!= '-'
    seq[length(primer)-findfirst(gapbool)+buffer+2:end]
end
function primer_trim(seq,phreds,primer; buffer = 3)
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq[1:length(primer)+buffer])
    gapbool = reverse(collect(a1)) .!= '-'
    seq[length(primer)-findfirst(gapbool)+buffer+2:end],phreds[length(primer)-findfirst(gapbool)+buffer+2:end]
end

export primer_peek
function primer_peek(seqs::Vector{String}; l = 20, N = 30, keep = 20)
    #Assumptions: No ambigs, and primer sequence goes "barcode+primer"
    starts = [s[1:l] for s in seqs]
    ends = [s[end-l+1:end] for s in seqs]
    pot_primers = countmap(vcat(starts,reverse_complement.(ends)))
    topN = reverse(sort([(pot_primers[k],k) for k in keys(pot_primers)]))[1:N];
    for i in enumerate(topN)
        println(i)
    end
    return reverse_complement.(sort(reverse_complement.([i[2] for i in topN[1:keep]])))
end
