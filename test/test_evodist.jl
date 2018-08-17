@testset "evodist" begin
    @testset "estimate_distance" begin
        out = estimate_distance("GTCGATCGACTAGCT", "GACGCTCGACTAGGT")
        @test out == 0.23963672893825572
    end

end