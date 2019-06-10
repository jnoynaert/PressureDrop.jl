@testset "Pressure and temp wrapper" begin

segments = 100
MDs = range(0, 5000, length = segments) |> collect
incs = repeat([0], inner = segments)
TVDs = range(0, 5000, length = segments) |> collect

well = Wellbore(MDs, incs, TVDs, 2.441)

model = WellModel(wellbore = well, roughness = 0.0006,
                    temperature_method = "Shiu", geothermal_gradient = 1.0, BHT = 200,
                    pressurecorrelation = HagedornAndBrown, WHP = 350 - pressure_atmospheric, dp_est = 25,
                    q_o = 100, q_w = 500, GLR = 1200, APIoil = 35, sg_water = 1.1, sg_gas = 0.8)

pressures = pressure_and_temp!(model)
temps = model.temperatureprofile

@test length(pressures) == length(temps) == segments
@test pressures[1] == 350 - pressure_atmospheric
@test pressures[end] ≈ (1068 - pressure_atmospheric) atol = 1
@test temps[end] == 200
@test temps[1] ≈ 181 atol = 1

end #testset for pressure & temp wrapper

#TODO: redo this using read_valves (as another implicit check) and use as an end-to-end test:
#=
end to end GL test:
segments = 100
MDs = range(0, 5000, length = segments) |> collect
incs = repeat([0], inner = segments)
TVDs = range(0, 5000, length = segments) |> collect

well = Wellbore(MDs, incs, TVDs, 2.441)

valves = GasliftValves([1813,2375,2885,3395], [1005,990,975,960], [0.073,0.073,0.073,0.073], [16,16,16,16])

model = WellModel(wellbore = well, roughness = 0.0006, valves = valves,
                    temperature_method = "Shiu", geothermal_gradient = 1.0, BHT = 200,
                    pressurecorrelation = HagedornAndBrown, WHP = 350 - pressure_atmospheric, dp_est = 25,
                    q_o = 100, q_w = 500, GLR = 1200, APIoil = 35, sg_water = 1.1, sg_gas = 0.8, CHP = 1000, naturalGLR = 1)

tubing_pressures, casing_pressures, valvedata = gaslift_model!(model, find_injectionpoint = true, dp_min = 100)

plot_gaslift(model.wellbore, tubing_pressures, casing_pressures, model.temperatureprofile, valvedata, nothing)
=#
