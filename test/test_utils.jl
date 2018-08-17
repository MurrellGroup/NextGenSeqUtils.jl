@testset "utils" begin
    #=
    @testset "dist_matrix" begin
        input = [2, 2, 2, 2]
        dist = dist_matrix(input, input; dist_met = (x,y)->abs(x-y))
        @test dist == [[0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0],  
                      [0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0]]
    end
    =#

    @testset "degap" begin
        input = "GGG-G-A-CT"
        degapped = degap(input)
        @test degapped == "GGGGACT"
    end
    
    @testset "dash_count" begin
        input = "GGG-G-A-CT"
        count = dash_count(input)
        @test count == 3
    end
    
    @testset "single_gap" begin
        input = "GGG-G-A-CT"
        dash = single_gap(input, 4)
        @test dash == true
    end
    
    @testset "logsum" begin
        out = logsum(1, 2)
        @test out == 2.313261687518223
    end
end