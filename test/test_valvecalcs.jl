include("../src/valvecalculations.jl")

@testset "Dome pressures" begin
    #testing with R-value of zero to ignore the adjustment to PTRO
    @test domepressure_downhole(600-14.65, 0, 180) ≈ 750 atol = 5 #function takes psig but test values are psia
    @test domepressure_downhole(727-14.65, 0, 120) ≈ 820 atol = 5
end


@testset "Thornhill-Craver" begin
    @test ThornhillCraver_gaspassage_simplified(900, 1100, 140, 16) ≈ 1127 atol=1 #using psig inputs
    @test ThornhillCraver_gaspassage_simplified(1200, 1100, 140, 16) == 0

    seats = [7*4, 32, 5*8, 3*16]
    Hernandez_results = [3.15, 4.11, 6.36, 9.39] .* 1000
    results = [ThornhillCraver_gaspassage(420 - pressure_atmospheric, 850 - pressure_atmospheric, 150, s, 0.7) for s in seats]
    @test all(abs.(Hernandez_results .- results) .<= Hernandez_results .* 0.02) #2% tolerance to account for changing C_ds that aren't easily available (test cases use Winkler C_ds)

    @test ThornhillCraver_gaspassage(1200, 1100, 140, 16, 0.7) == 0
end


@testset "Valve table" begin

#EHU 256H example using Weatherford method

MDs = [0,1813, 2375, 2885, 3395]
TVDs = [0,1800, 2350, 2850, 3350]
incs = [0,0,0,0,0]
id = 2.441

well = Wellbore(MDs, incs, TVDs, id)
valves = GasliftValves([1813,2375,2885,3395], [1005,990,975,960], [0.073,0.073,0.073,0.073], [16,16,16,16])

tubing_pressures = [150,837,850,840,831]
casing_pressures = 1070 .+ [0,53,70,85,100]
temps = [135,145,148,151,153]

vdata, active_valve_row = valve_calcs(valves = valves, well = well, sg_gas = 0.72, tubing_pressures = tubing_pressures, casing_pressures = casing_pressures, tubing_temps = temps, casing_temps = temps)

valve_table(vdata, active_valve_row) #implicit test

results = vdata[1:4, [5,13,12,4]] #PSC, PVC, PVO, PSO

expected_results =
[1050. 1103 1124 1071;
 1023 1092 1111 1042;
 996  1080 1099 1015;
 968  1068 1087 987]

@test all(abs.(expected_results .- results) .< (expected_results .* 0.01)) #1% tolerance due to using TCFs versus PVT-based dome correction, as well as rounding errors

end #testset for valve table
