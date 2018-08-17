@testset "FAD" begin
    
    @testset "FAD" begin
        out1, out2 = FAD(vcat([reverse(nl43env) for i in 1:10], [nl43env for i in 1:10]))
        @test sort(out1) == sort([nl43env, reverse(nl43env)] && out2 == [10.0, 10.0])
    end

end