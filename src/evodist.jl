"""
    get_pis(sequences)
"""

function get_pis(sequences)
    c= countmap((collect(string(sequences...))))
    return (c['A'], c['C'], c['G'], c['T']) ./ sum([c['A'], c['C'], c['G'], c['T']]);
end

"""
    get_transition_mat(seq1, seq2; include_dash_seqs = true)
"""
function get_transition_mat(seq1, seq2; include_dash_seqs = true)
    a1, a2 = kmer_seeded_align(seq1, seq2)
    
    indices = ['A', 'C', 'G', 'T']
    matrix = zeros(4,4)
    
    for i in 1:length(a1)
        if a1[i] == '-' || a2[i] == '-'
            continue
        else
            matrix[findfirst(indices, a1[i]), findfirst(indices, a2[i])] += 1
        end
    end
    return matrix 
end

"""
    estimate_distance(seq1, seq2)
"""
function estimate_distance(seq1, seq2)
    matrix = get_transition_mat(seq1, seq2)
    pi_a, pi_c, pi_g, pi_t = get_pis([seq1, seq2])
    pi_y = pi_t + pi_c
    pi_r = pi_a + pi_g
    
    s_1 = (matrix[2,4] + matrix[4,2])/sum(matrix)
    s_2 = (matrix[1,3] + matrix[3,1])/sum(matrix)
    v = (matrix[1,2] + matrix[2,1] + matrix[1,4] + matrix[4,1] + matrix[2,3] + matrix[3,2] + matrix[3,4] + matrix[4,3])/sum(matrix)
    
    a_1 = - log(1 - (pi_y * s_1)/(2 * pi_t * pi_c) - (v/(2 * pi_y) ))
    a_2 = - log(1 - (pi_r * s_2)/(2 * pi_a * pi_g) - (v/(2 * pi_r) ))
    b = - log(1 - (v/(2 * pi_y * pi_r)))
    
    d = ((2 * pi_t * pi_c)/pi_y)*(a_1 - (pi_r * b)) + ((2 * pi_a * pi_g)/pi_r)*(a_2 - (pi_y * b)) + (2 * pi_y * pi_r * b)
    
    return d
end
