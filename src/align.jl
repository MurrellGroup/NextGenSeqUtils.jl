"""
    nw_align(s1::String, s2::String; edge_reduction = 0.99)

Returns aligned strings using the Needleman-Wunch Algorithm (quadratic),
with end gaps penalized slightly less. edge_reduction is a multiplier (usually
less than one) on gaps on end of strings.
"""
function nw_align(s1::String, s2::String; edge_reduction = 0.99, mismatch_cost = -1.0)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    match_cost = 1.0

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
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = (i == length(s1arr)+1) ? edge_reduction : 1
            insMult = (j == length(s2arr)+1) ? edge_reduction : 1

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = findmax(scores)[2]
            arr[i, j] = scores[best]
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]

    # the trace endpoing will need to be generalized if we want to
    # allow overhang.
    trI, trJ = length(s1arr), length(s2arr)

    # First compute the trace running backwards. Initialized to the
    # maximum size. With unit penalties, is it possible to tell how
    # big this should be in advance?
    backtrace = Array{Int,1}(undef, length(s1arr) + length(s2arr))
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
    ali1arr = Array{Char,1}(undef, length(backtrace))
    ali2arr = Array{Char,1}(undef, length(backtrace))
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
    return join(ali1arr), join(ali2arr)
end


#--------Banded alignment---------

"""
    nw_align(s1::String, s2::String, banded::Float64)

Wrapper for `nw_align` and `banded_nw_align`. A larger `banded` value makes alignment slower
but more accurate.
"""
function nw_align(s1::String, s2::String, banded::Float64)
    if banded <= 0
        return nw_align(s1, s2)
    else
        return banded_nw_align(s1, s2, band_coeff = banded)
    end
end

# methods for converting from rectangular coordinates to diagonal band coordinates and back

"""Sets value in band at `i`, `j` to `val`, where `i` and `j` are in square matrix coords, `dim_diff` = ncols - nrows"""
function add_to_band!(band, val, i::Int, j::Int, bandwidth::Int, dim_diff::Int)
    # smallest coord is slice, largest is +/- offset from middle of slice + offset for diff in string lengths
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    band[ii, jj] = val
end

"""Returns value from band where `i` and `j` are in square matrix coords, `dim_diff` = ncols - nrows"""
function get_band_val(band, i::Int, j::Int, bandwidth::Int, dim_diff::Int)
    # smallest coord is slice, largest is +/- offset from middle of slice + offset for diff in string lengths
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    return band[ii, jj]
end

"""Checks if given square matrix coords are within width of band (ignore length)"""
function in_band(i, j, bandwidth, dim_diff)
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    return !( ii < 1 || jj < 1 || jj > 2*bandwidth + 1 + abs(dim_diff) )
end

"""
    banded_nw_align(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)

Like nw_align, but sub quadratic by only computing values within a band around the center diagonal.
One 'band' of radius 3 = (4,1), (3,1), (2,1), (1,1), (1,2), (1,3), (1,4), aka upside-down L shape.
band_coeff = 1 is sufficient to get same alignments as nw_align for 10% diverged sequences ~97% of the time;
increase this value for more conservative alignment with longer computation time.
Radius of band = `bandwidth` = `band_coeff` * sqrt(avg seq length)
"""
function banded_nw_align(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)
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
                    (get_band_val(arr, i-1, j, bandwidth, dim_diff)+(ins_cost*insMult)) : -Inf
            del = in_band(i, j-1, bandwidth, dim_diff) ?
                    (get_band_val(arr, i, j-1, bandwidth, dim_diff)+(del_cost*delMult)) : -Inf

            scores = [diag, del, ins]
            best = findmax(scores)[2]
            add_to_band!(arr, scores[best], i, j, bandwidth, dim_diff)
            add_to_band!(traceArr, best, i-1, j-1, bandwidth, dim_diff)
        end
    end
    alignedScore = get_band_val(arr, length(s1arr)+1, length(s2arr)+1, bandwidth, dim_diff)

    trI, trJ = length(s1arr), length(s2arr)
    backtrace = Array{Int,1}(undef, length(s1arr) + length(s2arr))
    btInd = 1
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = get_band_val(traceArr, trI, trJ, bandwidth, dim_diff)
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
        elseif backtrace[btInd] == 2
            trJ += -1
        elseif backtrace[btInd] == 3
            trI += -1
        else
            error("Bad trace value: $(backtrace[btInd]): ($trI, $trJ)")
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
    ali1arr = Array{Char,1}(undef, length(backtrace))
    ali2arr = Array{Char,1}(undef, length(backtrace))
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
    return join(ali1arr), join(ali2arr)
end


#----Triplet/Codon Alignment----

"""
    triplet_nw_align(s1::String, s2::String; edge_reduction = 0.99, boundary_mult = 2)

Returns alignment of two sequences where `s1` is a reference with reading frame to be preserved and `s2` is a query sequence.
`boundary_mult` adjusts penalties for gaps preserving the reading frame of `s1`.
This usually works best on range 0 to 3, higher values for more strongly enforced gaps aligned on
reading frame (divisible-by-3 indices)
"""
function triplet_nw_align(s1::String, s2::String; edge_reduction = 0.99, boundary_mult = 2)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical

    triple_ins_cost = -0.5
    triple_del_cost = -0.5
    triple_ins_score(qv) = triple_ins_cost
    triple_del_score(qv) = triple_del_cost

    mismatch_cost = -0.0
    match_cost = 1.0

    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1, :] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = (i == length(s1arr)+1) ? edge_reduction : 1
            insMult = (j == length(s2arr)+1) ? edge_reduction : 1
            del_bndry_mult = ((i-1) % 3 == 0) ? -boundary_mult : 1
            ins_bndry_mult = ((j-1) % 3 == 0) ? -boundary_mult : 1

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            # TODO: factor in qv scores
            s1qv = NaN
            s2qv = NaN
            triple_ins = i > 3 ? (arr[i-3, j]+(triple_ins_score(s1qv) * ins_bndry_mult)) : -Inf
            triple_del = j > 3 ? (arr[i, j-3]+(triple_del_score(s2qv) * del_bndry_mult)) : -Inf
            scores = [diag, del, ins, triple_del, triple_ins]
            best = findmax(scores)[2]
            arr[i, j] = scores[best]
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]

    trI, trJ = length(s1arr), length(s2arr)
    # First compute the trace running backwards
    backtrace = Array{Int,1}(undef, length(s1arr) + length(s2arr))
    btInd = 1
    len = 0
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = traceArr[trI, trJ]
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
            len += 1
        elseif backtrace[btInd] == 2
            trJ += -1
            len += 1
        elseif backtrace[btInd] == 3
            trI += -1
            len += 1
        elseif backtrace[btInd] == 4
            trJ += -3
            len += 3
        elseif backtrace[btInd] == 5
            trI += -3
            len += 3
        else
            error("uh oh: $(traceArr[trI, trJ])")
        end
        btInd += 1
    end
    # If you hit the boundaries not at the top left corner.
    while trI > 0
        backtrace[btInd] = 3
        btInd += 1
        trI += -1
        len += 1
    end
    while trJ > 0
        backtrace[btInd] = 2
        btInd += 1
        trJ += -1
        len += 1
    end

    backtrace = backtrace[1:btInd-1]
    ali1arr = Array{Char,1}(undef, len+3)
    ali2arr = Array{Char,1}(undef, len+3)
    ind1 = 1
    ind2 = 1
    i = 1
    for tr in reverse(backtrace)
        if tr == 1
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = s2arr[ind2]
            ind1 += 1
            ind2 += 1
            i += 1
        elseif tr == 2
            ali2arr[i] = s2arr[ind2]
            ali1arr[i] = '-'
            ind2 += 1
            i += 1
        elseif tr == 3
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = '-'
            ind1 += 1
            i += 1
        elseif tr == 4
            ali2arr[i:i+2] = s2arr[ind2:ind2+2]
            ali1arr[i:i+2] = ['-', '-', '-']
            ind2 += 3
            i += 3
        elseif tr == 5
            ali1arr[i:i+2] = s1arr[ind1:ind1+2]
            ali2arr[i:i+2] = ['-', '-', '-']
            ind1 += 3
            i += 3
        end
    end
    # to fix: weird boundary checking -- sometimes get off by 1 in length
    return join(ali1arr[1:i-1]), join(ali2arr[1:i-1])
end

codon_nw_align = triplet_nw_align


#-------Local Alignment--------

"""
    local_align(ref::String, query::String; mismatch_score = -1,
                match_score = 1, gap_penalty = -1,
                rightaligned=true, refend = false)

Aligns a query sequence locally to a reference. If true, `rightaligned` keeps the
right ends of each sequence in final alignment- otherwise they are trimmed;
`refend` keeps the beginning/left end of `ref`.
If you want to keep both ends of both strings, use nw_align.
For best alignments use the default score values.
"""
function local_align(ref::String, query::String; mismatch_score = -1,
                     match_score = 1, gap_penalty = -1,
                     rightaligned=true, refend = false)
    s1 = ref
    s2 = query
    # Populate Scores
    scoremat = zeros(length(s1)+1, length(s2)+1)
    scorematdirs = zeros(length(s1), length(s2))
    hiscore = 0
    hiscore_ij = (0, 0)
    for i in 2:size(scoremat)[1]
        for j in 2:size(scoremat)[2]
            # indexes indications: 1 = left and up, 2 = up, 3 = left, 4 = terminate (0)
            score_ij = [
                        scoremat[i-1, j-1] + (s1[i-1] == s2[j-1] ? match_score : mismatch_score),
                        scoremat[i-1, j] + gap_penalty,
                        scoremat[i, j-1] + gap_penalty,
                        0
                    ]
            # record index of max in matrix, track current value and location of high score
            ind = findmax(score_ij)[2]
            scoremat[i, j] = score_ij[ind]
            scorematdirs[i-1, j-1] = ind
            if score_ij[ind] > hiscore
                hiscore = score_ij[ind]
                hiscore_ij = (i-1, j-1)
            end
        end
    end

    # Traceback Path
    i, j = size(scorematdirs)
    if !rightaligned
        i, j = hiscore_ij
    end
    backtrace = Array{UInt8,1}(undef, length(s1) + length(s2))
    step = 0
    while i > 0 && j > 0
        direction = scorematdirs[i, j]
        if direction == 4
            break
        end
        step += 1
        backtrace[step] = direction
        if direction == 1
            # both chars
            i -= 1
            j -= 1
        elseif direction == 2
            # s1 char, s2 gap
            i -= 1
        else  # direction == 3
            # s1 gap, s2 char
            j -= 1
        end
    end
    # run along edges to end
    while refend && i > 0
        step += 1
        backtrace[step] = 2
        i -= 1
    end

    # Construct Strings
    ali1arr = Array{Char,1}(undef, step)
    ali2arr = Array{Char,1}(undef, step)
    ind1 = 0
    ind2 = 0
    for k in 1:step
        dir = backtrace[k]
        if dir == 1
            ali1arr[k] = s1[end-ind1]
            ali2arr[k] = s2[end-ind2]
            ind1 += 1
            ind2 += 1
        elseif dir == 2
            ali1arr[k] = s1[end-ind1]
            ali2arr[k] = '-'
            ind1 += 1
        elseif dir == 3
            ali2arr[k] = s2[end-ind2]
            ali1arr[k] = '-'
            ind2 += 1
        else
            error("error in local_align: alignment score matrix has invalid values")
        end
    end
    return join(reverse(ali1arr)), join(reverse(ali2arr))
end

#--------Kmer_seeded_align internals---------

"""
Sets value of key in dictionary to given index, unless the key already exists,
in which case value is set to -1.
"""
function unique_key(dicto::Dict{String, Int}, keyo::String, indo::Int)
    if haskey(dicto, keyo)
        dicto[keyo] = -1
    else
        dicto[keyo] = indo
    end
end

# Obsolete
function clean_unique_key(dicto::Dict{String, Int})
    for i in keys(dicto)
        if dicto[i] == -1
            delete!(dicto, i)
        end
    end
end

"""
Returns a list of indices of matches of unique words between `s1` and `s2`,
sorted by index of `s1`.
"""
function sorted_matches(s1, s2, wordlength, skip, aligncodons)
    word_dict1 = Dict{String, Int}()
    word_dict2 = Dict{String, Int}()
    # we can space one of these out, but not both
    for i in 1:skip:(length(s1)-(wordlength-1))  # bounds checked
        if aligncodons
            unique_key(word_dict1, translate_to_aa(s1[i:i+(wordlength-1)]), i)
        else
            unique_key(word_dict1, s1[i:i+(wordlength-1)], i)  # bounds checked
        end
    end
    for i in 1:(length(s2)-(wordlength-1))  # bounds checked
        if aligncodons
            unique_key(word_dict2, translate_to_aa(s2[i:i+(wordlength-1)]), i)
        else
            unique_key(word_dict2, s2[i:i+(wordlength-1)], i)  # bounds checked
        end
    end

    intersection = Dict{String, UInt8}()
    for (word, ind) in word_dict1
        if ind > 0 && haskey(word_dict2, word) && word_dict2[word] > 0
            intersection[word] = 0
        end
    end
    common = collect(keys(intersection))
    matches = zeros(Int, length(common), 2)
    for i in 1:length(common)
        matches[i, 1] = word_dict1[common[i]]
        matches[i, 2] = word_dict2[common[i]]
    end
    if length(matches)==0
        return matches
    else
        return sortslices(matches, by=x->(x[1]), dims=1)
    end
end

"""
Returns a list of indices of matches of unique words between `s1` and `s2`,
sorted by index of `s1`. Words match by amino acid encoding, in any reading frame.
"""
function sorted_aa_matches(str1, str2, wordlength)
    word_dict1 = Dict{String, Int}()
    word_dict2 = Dict{String, Int}()
    # get amino acid sequence in each reference frame
    s1s = generate_aa_seqs(str1)
    s2s = generate_aa_seqs(str2)
    for (offset, s1) in enumerate(s1s)
        for i in 1:(length(s1) - (div(wordlength, 3) - 1))
            # offset acounts for the shift in index due to being in a different reference frame.
            unique_key(word_dict1, s1[i:i+(div(wordlength, 3) - 1)], i*3 + offset - 3)
        end
    end
    for (offset, s2) in enumerate(s2s)
        for i in 1:(length(s2) - (div(wordlength, 3) - 1))
            unique_key(word_dict2, s2[i:i+(div(wordlength, 3) - 1)], i*3 + offset - 3)
        end
    end
    intersection = Dict{String, UInt8}()
    for (word, ind) in word_dict1
        if ind > 0 && haskey(word_dict2, word) && word_dict2[word] > 0
            intersection[word] = 0
        end
    end
    common = collect(keys(intersection))
    matches = zeros(Int, length(common), 2)
    for i in 1:length(common)
        matches[i, 1] = word_dict1[common[i]]
        matches[i, 2] = word_dict2[common[i]]
    end
    if length(matches)==0
        return matches
    else
        return sortslices(matches, by=x->(x[1]), dims=1)
    end
end

"""
Returns true if indices of word matches in second sequence (second column of
`matches`) is strictly increasing, else returns false.
"""
function matches_are_inconsistent(matches)
    # `matches` must already be sorted
    if size(matches)[1] <= 1
        return false
    end
    return minimum(matches[2:length(matches[:, 2]), 2] -
                   matches[1:(length(matches[:, 2]) - 1), 2]) < 1
end


"""
Some edge cases commonly arise where a kmer match starts before
the previous kmer match ends, but the two sequences still
mismatch.

To handle these, we make sure that either the word match in both sequences is spaced
`skip` after the previous word match, or both matches have spacing greater than
`wordlength`.
"""
function clean_matches(matches, wordlength, skip)
    # `matches` must already by sorted
    cleanedmatches = Vector{Int}[]
    push!(cleanedmatches, matches[1,:])
    currentvec = matches[1,:]
    for i in 2:length(matches[:, 2])
        # The condition below involves a kludge. It should really be
        # considering the case that skip == diff, rather than skip <
        # diff. Fix later.
        if !(skip < matches[i, 1] - currentvec[1] < wordlength ||
             (matches[i, 2] - currentvec[2] != skip &&
              matches[i, 2] - currentvec[2] < wordlength))
            push!(cleanedmatches, matches[i,:])
            currentvec = matches[i,:]
        end
    end
    return hcat(cleanedmatches...)'
end

"""
Returns ranges of overlapping word matches.
"""
function merge_overlapping(matches, wordlength, skip)
    range_inds = Vector{Int}[]
    push!(range_inds,[1, 1])
    current_range_ind_ind = 1
    for i in 2:length(matches[:, 2])
        # TODO: THESE BOUNDS (especially +wordlength) NEED TO BE CHECKED.
        # if ((matches[i, 1] < matches[i-1, 1] + (wordlength - 1) &&
        #      (matches[i, 2] < matches[i-1, 2] + (wordlength - 1))))
        if (matches[i, 1] - matches[i-1, 1] == skip) &&
            (matches[i, 2] - matches[i-1, 2] == skip)
            range_inds[current_range_ind_ind][2] = i
        else
            push!(range_inds,[i, i])
            current_range_ind_ind += 1
        end
    end
    return range_inds
end

"""
Returns actual word matches from indices of the matches.
"""
function get_matches(s1, s2, clean, range_inds, wordlength)
    s1matches = [clean[range_inds[i], 1] + [0, wordlength - 1]
                 for i in 1:length(range_inds)]
    s2matches = [clean[range_inds[i], 2] + [0, wordlength - 1]
                 for i in 1:length(range_inds)]
    match1 = [s1[i[1]:i[2]] for i in s1matches]
    match2 = [s2[i[1]:i[2]] for i in s2matches]
    return match1, match2
end

"""
Returns actual word mismatches (intervals between word matches) from indices of matches.
"""
function get_mismatches(s1, s2, matches, range_inds, wordlength)
    s1mismatches = [matches[(range_inds[i-1][2]):range_inds[i][1], 1] .+ [wordlength, -1]
                    for i in 2:length(range_inds)]
    s2mismatches = [matches[(range_inds[i-1][2]):range_inds[i][1], 2] .+ [wordlength,-1]
                    for i in 2:length(range_inds)]
    s1mismatches = vcat([[1, matches[1, 1]-1]],
                        s1mismatches,
                        [[matches[range_inds[length(range_inds)][2], 1]+wordlength, length(s1)]])
    s2mismatches = vcat([[1, matches[1, 2]-1]],
                        s2mismatches,
                        [[matches[range_inds[length(range_inds)][2], 2]+wordlength, length(s2)]])

    mismatch1 = [s1[i[1]:i[2]] for i in s1mismatches]
    mismatch2 = [s2[i[1]:i[2]] for i in s2mismatches]
    return mismatch1, mismatch2
end

"""
Find longest increasing subsequence of second column of given array.
Used to resolve bad orders of word matches while preserving as many matches as possible.
Because the first column (match locations in the first sequence) is sorted, corresponding
matches in the second column must be sorted, so we get the maximum number of such matches.
"""
function longest_incr_subseq(arr::Array{Int, 2})
    len = size(arr)[1]
    prevs = ones(Int, len)
    curr = ones(Int, len)
    maxlen = 0
    for i in 1:len
        # binary search for longest seq so far that's less than current value.
        lo = 1
        hi = maxlen
        while lo <= hi
            mid = Int(ceil((lo+hi)/2))
            if arr[curr[mid], 2] < arr[i, 2]
                lo = mid + 1
            else
                hi = mid - 1
            end
        end
        # update sequences seen and corresponding previous indices
        if lo > 1
            prevs[i] = curr[lo-1]
        else
            prevs[i] = 1
        end
        curr[lo] = i
        if lo > maxlen
            maxlen = lo
        end
    end
    # build subsequence from indices
    subseq = zeros(Int, maxlen, 2)
    k = curr[maxlen]
    for i in maxlen:-1:1
        subseq[i, :] = arr[k, :]
        k = prevs[k]
    end
    return subseq
end


#--------end of Kmer_seeded_align internals---------



"""
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
"""
function kmer_seeded_align(s1::String, s2::String;
                           wordlength = 30,
                           skip = 10,
                           aligncodons = false,
                           banded = 1.0,
                           debug::Bool = false)

    # ToDo (i think we tried this and it had insignificant speed up):
    # 1: Make this recurse. So instead fo calling nw_align on the set
    # of mismatches strings, it calls itself, but with a lower word
    # length. This will help for really noisy sequences.

    if (length(s1) < wordlength || length(s2) < wordlength || s1 == "" || s2 == "")
        return nw_align(s1, s2)
    end

    sorted = sorted_matches(s1, s2, wordlength, skip, aligncodons)
    if size(sorted)[1] == 0
        return nw_align(s1, s2, banded)
    end
    if matches_are_inconsistent(sorted)
        # try to extract an in-order subsequence to resolve inconsistency
        sorted = longest_incr_subseq(sorted)
        if matches_are_inconsistent(sorted)
            println("Notice: Word matching produced inconsistent ordering.",
                    " Consider a larger word size. Returning full DP alignment.",
                    "\nNote: this message shouldn't print with code changes. If it does,
                    check kner_seeded_align helper methods.")
            return nw_align(s1, s2, banded)
        end
    end

    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # get matches
    match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength)
    # get mismatches
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength)

    alignedStrings = collect(nw_align(mismatch1[1], mismatch2[1]))
    for i in 1:(length(match1))
        alignedStrings[1] = alignedStrings[1] * match1[i]
        alignedStrings[2] = alignedStrings[2] * match2[i]
        # al1, al2 = align(mismatch1[i+1], mismatch2[i+1],-1.0,-1.0, 20)
        al1, al2 = nw_align(mismatch1[i+1], mismatch2[i+1])
        alignedStrings[1] = alignedStrings[1] * al1
        alignedStrings[2] = alignedStrings[2] * al2
    end

    if debug
        if (degap(alignedStrings[1]) != degap(s1) || degap(alignedStrings[2]) != degap(s2))
            error("Aligned strings do not match original strings")
        end
    end
    return alignedStrings
end

"""
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
"""
function triplet_kmer_seeded_align(s1::String, s2::String;
                           wordlength = 30,
                           skip = 9,
                           boundary_mult = 2,
                           alignedcodons = true,
                           debug::Bool=false)
    if !(wordlength % 3 == skip % 3 == 0)
        error("wordlength and skip must be divisible by 3")
    end

    if s1 == "" || s2 == ""
        return triplet_nw_align(s1, s2, boundary_mult=boundary_mult)
    end
    sorted = sorted_matches(s1, s2, wordlength, skip, alignedcodons)
    if size(sorted)[1] == 0
        return triplet_nw_align(s1, s2, boundary_mult=boundary_mult)
    end
    if matches_are_inconsistent(sorted)
        println("Notice: Word matching produced inconsistent ordering.",
                " Consider a larger word size. Returning full DP alignment.")
        return triplet_nw_align(s1, s2, boundary_mult=boundary_mult)
    end

    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # get matches
    match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength)
    # get mismatches
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength)

    alignedStrings = collect(triplet_nw_align(mismatch1[1], mismatch2[1], boundary_mult=boundary_mult))
    for i in 1:(length(match1))
        alignedStrings[1] = alignedStrings[1] * match1[i]
        alignedStrings[2] = alignedStrings[2] * match2[i]
        # al1, al2 = align(mismatch1[i+1], mismatch2[i+1],-1.0,-1.0, 20)
        al1, al2 = triplet_nw_align(mismatch1[i+1], mismatch2[i+1], boundary_mult=boundary_mult)
        alignedStrings[1] = alignedStrings[1] * al1
        alignedStrings[2] = alignedStrings[2] * al2
    end

    debug = true
    if debug && (degap(alignedStrings[1]) != degap(s1) || degap(alignedStrings[2]) != degap(s2))
        error("Aligned strings do not match original strings")
    end
    return alignedStrings
end

"""
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
"""
function local_kmer_seeded_align(s1::String, s2::String;
                                   wordlength = 30,
                                   skip = 10,
                                   trimpadding = 100,
                                   debug::Bool=false)
    if s1 == "" || s2 == ""
        return local_align(s1, s2, rightaligned=false)
    end

    sorted = sorted_matches(s1, s2, wordlength, skip, false)
    if size(sorted)[1] == 0
        return local_align(s1, s2, rightaligned=false)
    end
    if matches_are_inconsistent(sorted)
        sorted = longest_incr_subseq(sorted)
        if matches_are_inconsistent(sorted)
            println("Notice: Word matching produced inconsistent ordering.",
                    " Consider a larger word size. Returning full DP alignment.")
            return local_align(s1, s2, rightaligned=false)
        end
    end

    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # get mismatches
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength)
    # partition into front, middle, end mismatches, where middle is between beginning and end matching words
    frontmismatch1, mismatch1, endmismatch1 = mismatch1[1], mismatch1[2:end-1], mismatch1[end]
    frontmismatch2, mismatch2, endmismatch2 = mismatch2[1], mismatch2[2:end-1], mismatch2[end]
    # get matches
    match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength)

    alignedStrings = ["", ""]
    for i in 1:(length(match1)-1)
        # append a match
        alignedStrings[1] = alignedStrings[1] * match1[i]
        alignedStrings[2] = alignedStrings[2] * match2[i]
        # append a mismatch
        al1, al2 = nw_align(mismatch1[i], mismatch2[i])
        alignedStrings[1] = alignedStrings[1] * al1
        alignedStrings[2] = alignedStrings[2] * al2
    end
    # append final match
    alignedStrings[1] = alignedStrings[1] * match1[end]
    alignedStrings[2] = alignedStrings[2] * match2[end]

    # trim query ends to trimpadding + length of ends of reference before local alignment to go fast
    if trimpadding > 0
        frontmismatch2 = frontmismatch2[max(1, end-length(frontmismatch1)-trimpadding):end]
        endmismatch2 = endmismatch2[1:end-length(endmismatch1)-trimpadding]
    end
    # do semiglobal alignment with beginning and ends and make them closer to the middle matches
    # reverse end mismatches to make local_align's default right alignment be left alignment
    fali1, fali2 = local_align(frontmismatch1, frontmismatch2, refend=true)
    eali1, eali2 = local_align(reverse(endmismatch1), reverse(endmismatch2), refend=true)
    eali1, eali2 = reverse(eali1), reverse(eali2)
    alignedStrings[1] = fali1*alignedStrings[1]*eali1
    alignedStrings[2] = fali2*alignedStrings[2]*eali2

    if debug
        if (degap(alignedStrings[1]) != degap(s1) || degap(alignedStrings[2]) != degap(s2))
            error("Aligned strings do not match original strings")
        end
    end
    return alignedStrings
end

loc_kmer_seeded_align = local_kmer_seeded_align

"""
    kmer_seeded_edit_dist(s1::String , s2::String;
                          wordlength = 30,
                          skip = 5,
                          aa_matches = false)

Computes levenshtein edit distance with speedups from only computing the dp scoring matrix between word matches.
If aa_matches = true, will attempt to find amino acid matches in any reference frame,
and add the nucleotide Hamming distance of these matches to Levenshtein distances of mismatches.
`skip` gives a necessary gap between searched-for words in `s1`.
For best results, use the default `wordlength` and `skip` values.
"""
function kmer_seeded_edit_dist(s1::String , s2::String;
                               wordlength = 30,
                               skip = 5,
                               aa_matches = false)
    if (length(s1) < wordlength || length(s2) < wordlength || s1 == "" || s2 == "")
        return levenshtein(s1, s2)
    end
    if aa_matches
        if wordlength % 3 != 0
            error("If aa_matches == true, wordlength must be divisible by 3")
        end
        sorted = sorted_aa_matches(s1, s2, wordlength)
        skip = 1
    else
        sorted = sorted_matches(s1, s2, wordlength, skip, false)
    end
    if size(sorted)[1] == 0
        return levenshtein(s1, s2)
    end
    if matches_are_inconsistent(sorted)
        # try to extract an in-order subsequence to resolve inconsistency
        sorted = longest_incr_subseq(sorted)
        if matches_are_inconsistent(sorted)
            println("Notice: Word matching produced inconsistent ordering.",
                    " Consider a larger word size. Returning full edit distance.")
            return levenshtein(s1, s2)
        end
    end
    trim = 0
    if aa_matches
        # trim word matches to avoid double counting some errors on ends of matches
        trim = max(div(wordlength, 5), 1)
        for i in 1:size(sorted)[1]
            sorted[i, 1] += trim
            sorted[i, 2] += trim
        end
    end
    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # 2*trim to accomodate trimmed word matches on each side of the word match
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength - 2*trim)
    matched_diffs = 0
    if aa_matches
        match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength - 2*trim)
        for i in 1:length(match1)
            for j in 1:length(match1[i])
                if match1[i][j] != match2[i][j]
                    matched_diffs += 1
                end
            end
        end
    end
    return matched_diffs + sum([levenshtein(mismatch1[i], mismatch2[i]) for i in 1:length(mismatch1)])
end


resolve_alignments(alignments::Array{String, 1}; mode = 1) = resolve_alignments(alignments[1], alignments[2], mode = mode)

"""
    resolve_alignments(ref::String, query::String; mode = 1)

Called on aligned strings. Resolves `query` with respect to `ref`.
`mode` = 1 for resolving single indels, `mode` = 2 for resolving single indels and codon insertions in query.
"""
function resolve_alignments(ref::String, query::String; mode = 1)
    if length(ref) != length(query)
        error("Arguments must be aligned")
    end
    for i in 1:length(ref)
        if ref[i] == query[i]
            continue
        end
        # 1%3 deletions, with gap of length 4 or something
#        if ((i%3 == 1 && i < length(query) && query[i+1] != '-') ||
#            (i%3 == 0 && query[i-1] != '-')) && single_mod_three_gap(query, i)
#            query = query[1:i-1] * "N" * query[i+1:end]
        if query[i] == '-' && !triple_gap(query, i)
            query = query[1:i-1] * "N" * query[i+1:end]
        # single insertions in query
        elseif single_gap(ref, i)
            query = query[1:i-1] * "-" * query[i+1:end]
        # single deletions in query
#        elseif single_gap(query, i)
#            query = query[1:i-1] * "N" * query[i+1:end]
        # copy codon of reference
        elseif mode == 2 && i <= length(ref)-3 && ref[i] == ref[i+1] == ref[i+2] == '-'
            query = query[1:i-1] * "---" * query[i+3:end]
        end
    end
    return degap(ref), degap(query)
end

"""
    align_reading_frames(clusters; k = 6, thresh = 0.03, verbose = false)

Takes `clusters` = [consensus_sequences, cluster_sizes], chooses references
out of consensuses that do not have stop codons in the middle, and makes all
consensus sequence reading frames agree. Returns resolved consensus seqs
(`goods`) along with filtered out consensus seqs that are >`thresh` divergent from
nearest reference (`bads`).
`k` = kmer size for computing kmer vectors of sequences.
"""
function align_reading_frames(clusters; k = 6, thresh = 0.03, verbose = false)
    consensuses = clusters[1]
    refs = String[]
    refkmers = []
    # find references
    for sq in consensuses
        tr = translate_to_aa(sq[1:end - length(sq)%3])
        if !('*' in tr[1:end-1])
            push!(refs, sq)
            push!(refkmers, kmer_count(sq, k))
        end
    end
    if length(refs) == 0
        return clusters, [[], []]
    end
    if verbose
        println("Using $(length(refs)) references to align reading frames")
    end
    goods = ([], [])
    bads = ([], [])
    # find closest ref and align
    for (i, sq) in enumerate(consensuses)
        sqkmer = kmer_count(sq, k)
        dists = zeros(length(refs))
        for (j, rfk) in enumerate(refkmers)
            dists[j] = corrected_kmer_dist(sqkmer, rfk, k=k)
        end
        minref_ind = indmin(dists)
        if dists[minref_ind] > thresh
            push!(bads[1], sq)
            push!(bads[2], clusters[2][i])
        else
            ali = triplet_kmer_seeded_align(refs[minref_ind], sq)
            #ali = triplet_kmer_seeded_align(refs[1], sq)
            _, newsq = resolve_alignments(ali...)
            # new reference frame may be really bad
            if (kmer_seeded_edit_dist(newsq, refs[minref_ind], aa_matches=true) > thresh*length(newsq))
                push!(bads[1], sq)
                push!(bads[2], clusters[2][i])
            else
                push!(goods[1], newsq)
                push!(goods[2], clusters[2][i])
            end
        end
    end
    if verbose
        println("Reading frame alignment: $(length(goods[1])) good sequences, $(length(bads[1])) bad sequences")
    end
    return goods, bads
end

align_reference_frames = align_reading_frames


"""
    local_edit_dist(s1::String, s2::String)

Returns the edit distance between two sequences after local alignment
"""
function local_edit_dist(s1::String, s2::String)
    str1, str2 = s1, s2
    if length(s1) > length(s2)
        str1, str2 = s2, s1
    end
    a, b = loc_kmer_seeded_align(str1, str2)
    dst = 0
    for i in 1:length(a)
        if a[i] != b[i] && a[i] != 'N' && b[i] != 'N'
            dst += 1
        end
    end
    return dst
end


"""
    affine_nw_align(s1::String, s2::String; gap_open = -2.0, gap_extend = -0.2, match_cost = 1.0, mismatch_cost = -1.0)
Performs full slow affine gap alignment. No tricks. No optimization. Use sparingly. Should likely be combined with kmer seeding.
"""
function NextGenSeqUtils.affine_nw_align(s1::String, s2::String;
    gap_open = -2.0,
    gap_extend = -0.2,
    match_cost = 1.0,
    mismatch_cost = -1.0,
    edge_reduction = 1.0)
   
    s1arr = collect(s1)  # vertical
    s1len = length(s1arr)
    s2arr = collect(s2)  # horizontal
    s2len = length(s2arr)
  
    M = zeros(s1len+1, s2len+1)
    IX = zeros(s1len+1, s2len+1)
    IY = zeros(s1len+1, s2len+1)

    IX[:,1] .= gap_open .+ edge_reduction * gap_extend * (0:s1len)
    IX[1,:] .= -Inf * [1:s2len+1;]
    IY[1,:] .= gap_open .+ edge_reduction * gap_extend * (0:s2len)
    IY[:,1] .= -Inf * [1:s1len+1;]
    M[1:end,1] .= gap_open .+ edge_reduction * gap_extend * (0:s1len) #Not sure about this.
    M[1,1:end] .= gap_open .+ edge_reduction * gap_extend * (0:s2len) #Not sure about this.
    
    traceM = zeros(Int, s1len+1, s2len+1)
    traceIX = zeros(Int, s1len+1, s2len+1)
    traceIY = zeros(Int, s1len+1, s2len+1)
    
    #ALSO UNSURE
    traceM[1,:] .+= 3 
    traceM[2:end,1] .+= 2
    
    traceIX[1,:] .+= 3 
    traceIY[2:end,1] .+= 2
    
    traceIY[1,:] .+= 3 
    traceIX[2:end,1] .+= 2
    
    for i in 2:s1len+1
        for j in 2:s2len+1
            diag_cost = 0.0
            if s1arr[i-1] == s2arr[j-1]
                diag_cost = match_cost
            else
                diag_cost = mismatch_cost
            end
            diagM = M[i-1, j-1] + diag_cost
            IX2M = IX[i-1, j-1] + diag_cost
            IY2M = IY[i-1, j-1] + diag_cost
            
            M[i,j],traceM[i,j] = findmax([diagM,IX2M,IY2M]) #1: stay in M. 2: come from IX. 3: come from IY
            
            M2IX = M[i-1, j] + gap_open
            IXextend = IX[i-1, j] + gap_extend * ifelse(i < s2len+1, 1.0, edge_reduction)
            IX[i,j],traceIX[i,j] = findmax([M2IX,IXextend]) #1: gap open. 2: gap extend
            
            M2IY = M[i, j-1] + gap_open
            IYextend = IY[i, j-1] + gap_extend * ifelse(i < s1len+1, 1.0, edge_reduction)
            IY[i,j],traceIY[i,j] = findmax([M2IY,-Inf,IYextend]) #1: gap open. 3: gap extend
              
        end
    end
    
    rev_arr1 = Char[]
    rev_arr2 = Char[]
    sizehint!(rev_arr1, s1len+s2len)
    sizehint!(rev_arr2, s1len+s2len)
    
    mats = [traceM,traceIX,traceIY]
    x_i = s1len+1
    y_i = s2len+1
    m_i = argmax([M[x_i,y_i],IX[x_i,y_i],IY[x_i,y_i]])
    while x_i-1 > 0 &&  y_i-1 > 0
        next_m_i = mats[m_i][x_i,y_i]
        if m_i == 1
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,s2arr[y_i-1])
            
            x_i -= 1
            y_i -= 1
        elseif m_i == 2
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,'-')
            x_i -= 1
        elseif m_i == 3
            push!(rev_arr1,'-')
            push!(rev_arr2,s2arr[y_i-1])
            y_i -= 1
        end
        m_i = next_m_i
    end
    while x_i-1 > 0
        push!(rev_arr1,s1arr[x_i-1])
        push!(rev_arr2,'-')
        x_i -= 1
    end
    while y_i-1 > 0
        push!(rev_arr1,'-')
        push!(rev_arr2,s2arr[y_i-1])
        y_i -= 1
    end
    return join(reverse(rev_arr1)),join(reverse(rev_arr2))
end

"""
A slow profile cost function. Compares two profile "columns" and scores how similar they are.
If you design you own, try and make them return 1 for perfectly similar and -1 for completely different.
"""
function profile_cost(p1::Array{Tuple{Char,Float64},1},p2::Array{Tuple{Char,Float64},1})
    cost = 0
    for i in p1
        for j in p2
            if i[1]==j[1]
                cost += i[2]*j[2]
            else
                cost = cost - i[2]*j[2]
            end
        end
    end
    return cost
end

"""
A pointless looking function that returns what a profile "gap" element should look like.
"""
function gap_elem(typo::Array{Array{Tuple{Char,Float64},1},1})
    return [('!',1.0)]
end

"""
    profile_affine_align(s1, s2, distfunc; gap_open = -2.0,gap_extend = -0.2)

Aligns two profiles, using an affine gap alignment strategy, inserting gap_elem(s1) gap elements wherever indicated.
"distfunc" should be "profile_cost" if you're using seqs2profile() to get profiles from sequence alignments.
ToDo: Add edge_reduction factor into profile_affine_align 

Usage:
nom1,seqs1 = read_fasta_with_names("alignment1.fasta");
nom2,seqs2 = read_fasta_with_names("alignment2.fasta");
prof1 = seqs2profile(seqs1);
prof2 = seqs2profile(seqs2);
alip1,alip2 = profile_affine_align(prof1,prof2,profile_cost);
"""
function profile_affine_align(s1, s2, distfunc;
    gap_open = -2.0,
    gap_extend = -0.2)
   
    s1arr = collect(s1)  # vertical
    s1len = length(s1arr)
    s2arr = collect(s2)  # horizontal
    s2len = length(s2arr)
  
    M = zeros(s1len+1, s2len+1)
    IX = zeros(s1len+1, s2len+1)
    IY = zeros(s1len+1, s2len+1)

    IX[:,1] .= gap_open .+ gap_extend * (0:s1len)
    IX[1,:] .= -Inf * [1:s2len+1;]
    IY[1,:] .= gap_open .+ gap_extend * (0:s2len)
    IY[:,1] .= -Inf * [1:s1len+1;]
    M[1:end,1] .= gap_open .+ gap_extend * (0:s1len) #Not sure about this.
    M[1,1:end] .= gap_open .+ gap_extend * (0:s2len) #Not sure about this.
    
    traceM = zeros(Int, s1len+1, s2len+1)
    traceIX = zeros(Int, s1len+1, s2len+1)
    traceIY = zeros(Int, s1len+1, s2len+1)
    
    #ALSO UNSURE
    traceM[1,:] .+= 3 
    traceM[2:end,1] .+= 2
    
    traceIX[1,:] .+= 3 
    traceIY[2:end,1] .+= 2
    
    traceIY[1,:] .+= 3 
    traceIX[2:end,1] .+= 2
    
    for i in 2:s1len+1
        for j in 2:s2len+1
            diag_cost = profile_cost(s1arr[i-1],s2arr[j-1])
            
            diagM = M[i-1, j-1] + diag_cost
            IX2M = IX[i-1, j-1] + diag_cost
            IY2M = IY[i-1, j-1] + diag_cost
            
            M[i,j],traceM[i,j] = findmax([diagM,IX2M,IY2M]) #1: stay in M. 2: come from IX. 3: come from IY
            
            M2IX = M[i-1, j] + gap_open
            IXextend = IX[i-1, j] + gap_extend
            IX[i,j],traceIX[i,j] = findmax([M2IX,IXextend]) #1: gap open. 2: gap extend
            
            M2IY = M[i, j-1] + gap_open
            IYextend = IY[i, j-1] + gap_extend
            IY[i,j],traceIY[i,j] = findmax([M2IY,-Inf,IYextend]) #1: gap open. 3: gap extend
              
        end
    end
    
    rev_arr1 = typeof(s1)()
    rev_arr2 = typeof(s1)()
    sizehint!(rev_arr1, s1len+s2len)
    sizehint!(rev_arr2, s1len+s2len)
    
    mats = [traceM,traceIX,traceIY]
    x_i = s1len+1
    y_i = s2len+1
    m_i = argmax([M[x_i,y_i],IX[x_i,y_i],IY[x_i,y_i]])
    while x_i-1 > 0 &&  y_i-1 > 0
        next_m_i = mats[m_i][x_i,y_i]
        if m_i == 1
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,s2arr[y_i-1])
            
            x_i -= 1
            y_i -= 1
        elseif m_i == 2
            push!(rev_arr1,s1arr[x_i-1])
            push!(rev_arr2,gap_elem(s1))
            x_i -= 1
        elseif m_i == 3
            push!(rev_arr1,gap_elem(s1))
            push!(rev_arr2,s2arr[y_i-1])
            y_i -= 1
        end
        m_i = next_m_i
    end
    while x_i-1 > 0
        push!(rev_arr1,s1arr[x_i-1])
        push!(rev_arr2,gap_elem(s1))
        x_i -= 1
    end
    while y_i-1 > 0
        push!(rev_arr1,gap_elem(s1))
        push!(rev_arr2,s2arr[y_i-1])
        y_i -= 1
    end
    return (reverse(rev_arr1)),(reverse(rev_arr2))
end

"""
    seqs2profile(seqs::Vector{String})

Converts a vector of aligned sequences into a "profile", where each column is represented by a vector of tuples.
The first tuple element is the character, and the second is the frequency.
"""
function seqs2profile(seqs::Vector{String})
    if length(union(length.(seqs))) != 1
        error("Sequences must be aligned (ie. be the same length)")
    end
    prof = Array{Array{Tuple{Char,Float64},1},1}([])
    sizehint!(prof,length(seqs[1]))
    for i in 1:length(seqs[1])
        char_dic = proportionmap([s[i] for s in seqs])
        push!(prof , [(k,char_dic[k]) for k in keys(char_dic)])
    end
    return prof
end


"""
    merge_alignments(seqs1::Vector{String},seqs2::Vector{String}; new_gap_char = '-')

Merges two sequence alignments using profile alignment. Not terribly fast, but useful.
Return type similar to input type, but with gaps inserted.

Usage example:
nom1,seqs1 = read_fasta_with_names("alignment1.fasta");
nom2,seqs2 = read_fasta_with_names("alignment2.fasta");
ali_s1,ali_s2 = merge_alignments(seqs1,seqs2);
write_fasta("merged.fasta",vcat(ali_s1,ali_s2),names = vcat(nom1,nom2));
"""
function merge_alignments(seqs1::Vector{String},seqs2::Vector{String}; new_gap_char = '-')
    prof1 = seqs2profile(seqs1);
    prof2 = seqs2profile(seqs2);
    
    alip1,alip2 = profile_affine_align(prof1,prof2,profile_cost);
    
    PROF_GAP = gap_elem(prof1)
    char_matrix = fill(new_gap_char,length(seqs1)+length(seqs2),length(alip1));

    coord1 = 1
    coord2 = 1
    len1 = length(seqs1)
    for s_col in 1:length(alip1)
        if alip1[s_col] != PROF_GAP
            char_matrix[1:len1,s_col] = [s[coord1] for s in seqs1]
            coord1 += 1
        end
        if alip2[s_col] != PROF_GAP
            char_matrix[len1+1:end,s_col] = [s[coord2] for s in seqs2]
            coord2 += 1
        end
    end
    return [join(char_matrix[i,:]) for i in 1:length(seqs1)],[join(char_matrix[i,:]) for i in length(seqs1)+1:length(seqs1)+length(seqs2)]
end
