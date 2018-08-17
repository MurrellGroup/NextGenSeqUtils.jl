@testset "kmers" begin
    @testset "kmer_count" begin
        input = "AAAAAC"
        bin = kmer_count(input, 5)
        @test bin[1] == 1
        @test bin[2] == 1
    end

    @testset "corrected_kmer_dist" begin
        input1 = [5.0, 5.0, 3.0, 2.0]
        input2 = [5.0, 6.0, 6.0, 2.0]
        k = 2
        ret = corrected_kmer_dist(input1, input2; k = 2)
        
        @test ret == sqeuclidean(input1, input2)/ (k*(sum(input1) + sum(input2)))
    end

end