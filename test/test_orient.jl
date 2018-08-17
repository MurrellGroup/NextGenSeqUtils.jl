@testset "orient" begin
    
    @testset "orient_strands" begin
        seqs = ["GTCGATCGACTAGCTGCATGACTGACATCGACATCGACGGAGCATGACTAGGACGACGAGCATC", 
                        "GTCGATCGACTAGCTGCATGACTGACATCGACATCGACGGAGCATGACTAGGACGACGAGCATC", 
                        "GTCGATCGACTAGCTGCATGACTGACATCGACATCGACGGAGCATGACTAGGACGACGAGCATCAAAAAAAAAAAAAAAAAAAAA"]  
        refs = "GTCGAT"
        names = ["s1", "s2", "s3"]
        phreds = [[Int8(30) for i in 1:length(seq)] for seq in seqs]
        oriented, newphreds, newnames = orient_strands(seqs, phreds, names, refs)
       
        @test oriented == ["GATGCTCGTCGTCCTAGTCATGCTCCGTCGATGTCGATGTCAGTCATGCAGCTAGTCGATCGAC",
                           "GATGCTCGTCGTCCTAGTCATGCTCCGTCGATGTCGATGTCAGTCATGCAGCTAGTCGATCGAC",
                           "TTTTTTTTTTTTTTTTTTTTTGATGCTCGTCGTCCTAGTCATGCTCCGTCGATGTCGATGTCAGTCATGCAGCTAGTCGATCGAC"]
    end

end