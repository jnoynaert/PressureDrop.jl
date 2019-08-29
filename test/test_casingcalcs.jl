include("../src/casingcalculations.jl")

@testset "Casing pressure segment drop" begin
    sg_gas = 0.7
    P_pc = HankinsonWithWichertPseudoCriticalPressure(sg_gas, 0, 0)
    _, T_pc, _ = HankinsonWithWichertPseudoCriticalTemp(sg_gas, 0, 0)
    ΔP = calculate_casing_pressuresegment(300, 10, (100+180)/2, 4500, sg_gas, KareemEtAlZFactor, P_pc, T_pc)
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
                                        CHP = 300 - pressure_atmospheric, sg_gas = 0.7, dp_est = 10)

    @test pressures[end] ≈ (332.7 - pressure_atmospheric) atol = 5

    #model = WellModel(wellbore = testwell, temperatureprofile = temps, CHP = 300 - pressure_atmospheric, sg_gas_inj = 0.7, dp_est_inj = 10)
    #pressures2 = casing_traverse_topdown(model)
    #@test pressures2[end] ≈ (332.7 - pressure_atmospheric) atol = 5
end
