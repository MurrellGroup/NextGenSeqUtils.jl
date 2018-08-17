using NextGenSeqUtils
using Distances

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "nextgensequtils" begin
    @testset "align" begin
        include("test_align.jl")
    end

    @testset "demux" begin
        include("test_demux.jl")
    end

    @testset "evo_dist" begin
        include("test_evodist.jl")
    end

    @testset "FAD" begin
        include("test_FAD.jl")
    end

    @testset "hmm" begin
        include("test_hmm.jl")
    end

    @testset "kmers" begin
        include("test_kmers.jl")
    end

    @testset "orient" begin
        include("test_orient.jl")
    end

    @testset "phreds" begin
        include("test_phreds.jl")
    end

    @testset "utils" begin
        include("test_utils.jl")
    end
    
end
