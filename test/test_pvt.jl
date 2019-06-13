include("../src/pvtproperties.jl")


@testset "Gas PVT" begin

    @test LeeGasViscosity(0.7, 2000.0, 160.0, 0.75) ≈ 0.017 atol = 0.001
    @test HankinsonWithWichertPseudoCriticalTemp(0.65, 0.05, 0.08) .- (370.2, 351.4, 18.09) |> x -> all(abs.(x) .<= 1)
    @test HankinsonWithWichertPseudoCriticalPressure(0.65, 0.05, 0.08) ≈ 634.0 atol = 2
    @test PapayZFactor(634, 351.4, 600.0, 100.0) ≈ 0.921 atol = 0.01
    @test KareemEtAlZFactor(663.29, 377.59, 2000, 150) ≈ 0.8242 atol = 0.05
    @test KareemEtAlZFactor_simplified(663.29, 377.59, 2000, 150) ≈ 0.8242 atol = 0.15
    @test gasVolumeFactor(1200, 0.85, 200) ≈ 0.013 atol = 0.001
    @test gasDensity_insitu(0.916, 0.883, 346.6, 80.3) ≈ 1.79 atol = 0.01

end


@testset "Oil PVT" begin

    @test StandingSolutionGOR(30.0, 0.6, 800.0, 120.0) ≈ 121.4 atol = 1
    @test StandingSolutionGOR(41.06, 0.916, 346.6, 80.3) ≈ 109.7 atol = 1
    @test StandingOilVolumeFactor(30.0, 0.6, 121.4, 800.0, 120.0) ≈ 1.07 atol = 0.05
    @test oilDensity_insitu(41.06,  0.916,  109.7,  1.05) ≈ 49.9 atol = 0.5
    @test BeggsAndRobinsonDeadOilViscosity( 37.9,  120) ≈ 4.05 atol = 0.05
    @test GlasoDeadOilViscosity(37.9, 120) ≈ 2.30 atol = 0.1
    @test ChewAndConnallySaturatedOilViscosity(2.3, 769.0) ≈ 0.638 atol = 0.1

end


@testset "Water PVT" begin

    @test waterDensity_stb(1.12) ≈ 69.9 atol = 0.1
    @test GouldWaterVolumeFactor(50.0, 120.0) ≈ 1.01 atol = 0.05

end

@testset "interfacial tension" begin

    #@test gas_oil_interfacialtension
    #@test gas_water_interfacialtension

end
