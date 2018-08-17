@testset "align" begin
    @testset "nw_align" begin
        out1, out2 = nw_align("GTCGATCGACTAGCTA","GTCGATGCGCCTGCTA")
        @test out1 == "GTCGAT-CGACTAGCTA"
        @test out2 == "GTCGATGCGCCT-GCTA"
    end

    @testset "triplet_nw_align" begin
        out1, out2 = triplet_nw_align("GTCGATCGACTAGCT","GTCGATAGTCGACTAGCT")
        @test out1 == "GTCGAT---CGACTAGCT"
        @test out2 == "GTCGATAGTCGACTAGCT"
    end
    
end