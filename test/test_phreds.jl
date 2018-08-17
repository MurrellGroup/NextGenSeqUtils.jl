@testset "phreds" begin

    @testset "phred_to_p" begin
        out = phred_to_p(Int8(20))
        @test out == 0.01
    end

    @testset "p_to_phred" begin
        out = p_to_phred(Float64(0.01))
        @test out == 20
    end

end