include("../src/casingcalculations.jl")

@testset "Casing pressure segment drop" begin
    ΔP = casing_pressuresegment_topdown(300, 10, (100+180)/2, 0, 4500, 0.7, 0, 0,
                                        HankinsonWithWichertPseudoCriticalPressure, HankinsonWithWichertPseudoCriticalTemp, KareemEtAlZFactor)
    @test ΔP ≈ 332.7-300 atol = 1
end


@testset "Casing pressure traverse" begin
    md = [0, 2250, 4500]
    tvd = md
    inc = [0, 0, 0]
    id = 2.441

    temps = [100., 140., 180.]

    testwell = Wellbore(md, inc, tvd, id)

    pressures = casing_traverse_topdown(wellbore = testwell, temperatureprofile = temps,
                                        CHP = 300, sg_gas = 0.7, dp_est = 10)

    @test pressures[end] ≈ 332.7 atol = 5
end
