@testset "demux" begin
    @testset "IUPAC_equals" begin
        out = IUPAC_equals(toIUPACnum("N")[1],toIUPACnum("A")[1])
        @test out == true
    end
end