@testset "Wellbore object" begin

md_bad = [-1.,2.,3.,4.]
md_good = [1.,2.,3.,4.]
inc = [1.,2.,3.,4.]
tvd_bad = [-1.,2.,3.,4.]
tvd_good = [1.,2.,3.,4.]
id = [1.,1.,1.,1.]

#implicit test for adding the leading 0,0 survey point:
w = Wellbore(md_good, inc, tvd_good, id)

@test w.id[1] == w.id[2]

#implicit test for allowing negatives:
Wellbore(md_bad, inc, tvd_bad, id, true)

try
    Wellbore(md_bad, inc, tvd_good, id)
catch e
    @test e isa Exception
end

try
    Wellbore(md_good, inc, tvd_bad, id)
catch e
    @test e isa Exception
end

end #testset for Wellbore object


@testset "Pressure and temp wrapper" begin

segments = 100
MDs = range(0, 5000, length = segments) |> collect
incs = repeat([0], inner = segments)
TVDs = range(0, 5000, length = segments) |> collect

well = Wellbore(MDs, incs, TVDs, 2.441)

pressures, temps = pressure_and_temp(well = well, roughness = 0.0006,
                                    temperature_method = "Shiu", geothermal_gradient = 1.0, BHT = 200,
                                    pressurecorrelation = HagedornAndBrown, WHP = 350, dp_est = 25,
                                    q_o = 100, q_w = 500, GLR = 1200, APIoil = 35, sg_water = 1.1, sg_gas = 0.8)

@test length(pressures) == length(temps) == segments
@test pressures[1] == 350
@test pressures[end] ≈ 1068 atol = 1
@test temps[end] == 200
@test temps[1] ≈ 181 atol = 1

end #testset for pressure & temp wrapper


@testset "Valve table" begin

#EHU 256H example using Weatherford method

MDs = [0,1813, 2375, 2885, 3395]
TVDs = [0,1800, 2350, 2850, 3350]
incs = [0,0,0,0,0]
id = 2.441

well = Wellbore(MDs, incs, TVDs, id)
valves = GasliftValves([1813,2375,2885,3395], [1005,990,975,960], [0.073,0.073,0.073,0.073], [16,16,16,16])

tubing_pressures = 14.7 .+ [150,837,850,840,831]
casing_pressures = 1070 .+ 14.7 .+ [0,53,70,85,100]
temps = [135,145,148,151,153]

vdata = valve_calcs(valves, well, 0.72, tubing_pressures, casing_pressures, temps, temps)

valve_table(vdata) #implicit test

results = vdata[1:4, [5,13,12,4]] #PSC, PVC, PVO, PSO

expected_results =
[1050. 1103 1124 1071;
 1023 1092 1111 1042;
 996  1080 1099 1015;
 968  1068 1087 987]

expected_results[:,2:3] = expected_results[:,2:3] .+ 14.7

@test all(abs.(expected_results .- results) .< (expected_results .* 0.01)) #1% tolerance due to using TCFs versus PVT-based dome correction, as well as rounding errors

end #testset for valve table
