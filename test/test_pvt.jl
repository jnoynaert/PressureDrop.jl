using Test

include("pvptproperties.jl")


@testset "Gas PVT" begin

    @test LeeGasViscosity(0.7, 2000.0, 160.0, 0.75) ≈ 0.017 atol = 0.001
    @test HankinsonWithWichertPseudoCriticalTemp(0.65, 0.05, 0.08) .- (370.2, 351.4, 18.09) |> x -> all(abs.(x) .<= 1)
    @test HankinsonWithWichertPseudoCriticalPressure(0.65, 0.05, 0.08) ≈ 634.0 atol = 2
    @test PapayZFactor(634, 351.4, 600.0, 100.0) ≈ 0.921 atol = 0.01
    @test gasVolumeFactor(1200, 0.85, 200) ≈ 0.013 atol = 0.001
    #TODO: @test gasDensity_insitu

end


@testset "Oil PVT" begin

    @test StandingSolutionGOR(30.0, 0.6, 800.0, 120.0) ≈ 121.4 atol = 1
    @test StandingOilVolumeFactor(30.0, 0.6, 121.4, 800.0, 120.0) ≈ 1.07 atol = 0.05
    #TODO: @test oilDensity_insitu(APIoil,  specificGravityGas,  solutionGOR,  oilVolumeFactor)
    @test BeggsAndRobinsonDeadOilViscosity( 37.9,  120) ≈ 4.05 atol = 0.05 #TODO: find another test result. I don't trust the source of this one.
    @test GlasoDeadOilViscosity(37.9, 120) ≈ 2.30 atol = 0.1
    @test ChewAndConnallySaturatedOilViscosity(2.3, 769.0) ≈ 0.638 atol = 0.1

end


@testset "Water PVT" begin

    @test waterDensity_stb(1.12) ≈ 69.9 atol = 0.1
    @test GouldWaterVolumeFactor(50.0, 120.0) ≈ 1.01 atol = 0.05
    @test assumedWaterViscosity == 1.0
    #TODO: @test waterDensity_insitu

end
