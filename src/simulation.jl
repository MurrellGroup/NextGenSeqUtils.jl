## Sequence simulation functions

#-------Simple simulations-------#

const CHAR_NUC_ALPH = ['A', 'C', 'G', 'T']
const NUCALPH = ["A", "C", "G", "T"]

const SIMPLE_ERRORS = vcat(NUCALPH, NUCALPH, NUCALPH, NUCALPH,
                           [join([NUCALPH[(i - 1) % 4 + 1],
                                  NUCALPH[Int64(floor((i + -1) / 4) % 4 + 1)]])
                            for i in 1:16],
                           ["" for i in 1:16])

"""
    flip(p, t, f)

Flip a biased coin.
"""
function flip(p, t, f)
    if rand() > p
        return t
 else return f
    end
end

"""
    simple_gen_seq(n::Int)

Generates a random uniform sequence of nucleotides.
"""
function simple_gen_seq(n::Int)
    return(join(sample(NUCALPH, n)))
end

"""
    simple_evolve(refseq, err_rate)

Evolves a sequence, uniformly.
"""
function simple_evolve(refseq, err_rate)
    return join([flip(err_rate, refseq[i], sample(SIMPLE_ERRORS, 1)[1])
                 for i in 1:length(refseq)])
end

"""
    fixed_diff_evolve(template::String, n_diffs::Int64)

Creates a fixed number of mutations of a sequence.
"""
function fixed_diff_evolve(template::String, n_diffs::Int64)
    diff_pos = sample(1:length(template), n_diffs, replace=false)
    char_arr = collect(template)
    for i in 1:n_diffs
        char_arr[diff_pos[i]] = sample([x for x in CHAR_NUC_ALPH if x != char_arr[diff_pos[i]]])
    end
    return join(char_arr)
end


#-------Functions for simulating amplicon sequences with a PacBio error model--------

"""
    run_length_encode(x::String)

Run length encoding of a string.
"""
run_length_encode(x::String) = StatsBase.rle([c for c in x])

"""
    run_length_decode(chars, lengths)

Run length decoding of a string.
"""
run_length_decode(chars, lengths) = String(StatsBase.inverse_rle(chars, lengths))

"""
    pb_error_inflation(old_length::Int64)

Homopolymer length-to-error-rate scaling based on PacBio sequencing.

This function may need to be tweaked to approximate how the error
rate increases with HP length.
"""
function pb_error_inflation(old_length::Int64)
    return old_length ^ 1.5
end

"""
    length_error_func(old_length::Int64; rate = 0.002)

Takes a true homopolymer length, and returns an observed homopolymer length based on PacBio sequencer error model.
"""
function length_error_func(old_length::Int64; rate = 0.002)
    direction = rand() < 0.5
    amount = rand(Poisson(rate * pb_error_inflation(old_length)))
    if direction
        return old_length + amount
    else
        return max(old_length - amount, 0)
    end
end

"""
    pb_seq_sim(template::String, rate::Float64; with_qvs = false)

Performs a sequence simulation from a template, specifying a target error rate, based on PacBio sequencing error model.
"""
function pb_seq_sim(template::String, rate::Float64; with_qvs = false)
    mut_prob = rate / 4
    chars, lengths = run_length_encode(template)
    # This 1.55 will need to be recalibrated if the HP distribution,
    # or the pb_error_inflation functions change.
    calib = 1.55
    new_lengths = [length_error_func(l, rate=rate/calib) for l in lengths]
    new_seq = run_length_decode(chars, new_lengths)
    new_string_seq = join([ifelse(rand() < mut_prob, CHAR_NUC_ALPH[rand(1:4)], i) for i in new_seq])

    if with_qvs
        # ADD IN ERROR PROB GENERATION. CONSIDER ADDING A
        # "DISCRIMINATION" PARAMETER THAT CONTROLS HOW INFORMATIVE THE
        # ERROR PROBS ARE OF ERROR LOCATIONS.
        tot_length = sum(new_lengths)
        baseErrorProbs = Array{Float64}(tot_length)
        current = 1
        for i in 1:length(new_lengths)
            for j in 1:new_lengths[i]
                baseErrorProbs[current] = (rate / 4) +
                    (rate / calib) *
                    (pb_error_inflation(new_lengths[i]) / new_lengths[i])
                current += 1
            end
        end

        perturbedLogs = log(baseErrorProbs) + randn(length(baseErrorProbs)) .* 0.1
        perturbed_qvs = e .^ perturbedLogs

        # This next line is roughly inspired by inverting poisson
        # parameter inference. Good luck.
        # totalErrorTweak = rand(Gamma(1+Int64(round(rate*length(perturbed_qvs))), 1))/length(perturbed_qvs)
        # perturbed_qvs = perturbed_qvs.*(totalErrorTweak/mean(perturbed_qvs))
        return [new_string_seq, perturbed_qvs]
    else
        return new_string_seq
    end
end

"""
    env_error_rates(n)

Draws from the error rate distribution typically seen in P5 envelope
sequence data.
"""
function env_error_rates(n)
    rand(Gamma(2, 0.0017), n)
end

"""
    env_pb_seq_sim(template::String, n::Int64; with_qvs = false)

Simulated PacBio reads from amplicons that have envelope-like error profiles.
"""
function env_pb_seq_sim(template::String, n::Int64; with_qvs = false)
    error_rates = env_error_rates(n)
    return [pb_seq_sim(template, err, with_qvs = with_qvs) for err in error_rates]
end
