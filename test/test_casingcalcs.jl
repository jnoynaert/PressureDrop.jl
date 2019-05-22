include("../src/casingcalculations.jl")

@testset "Casing pressure drop" begin
    ΔP = casing_pressure_segment_topdown(300, 10, (100+180)/2, 0, 4500, 0.7, 0, 0,
                                        HankinsonWithWichertPseudoCriticalPressure, HankinsonWithWichertPseudoCriticalTemp, KareemEtAlZFactor)
    @test ΔP ≈ 332.7-300 atol = 1
end
