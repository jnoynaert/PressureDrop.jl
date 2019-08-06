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


@testset "Gas lift wrappers" begin

segments = 100
MDs = range(0, 5000, length = segments) |> collect
incs = repeat([0], inner = segments)
TVDs = range(0, 5000, length = segments) |> collect

well = Wellbore(MDs, incs, TVDs, 2.441)

testpath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata")
valvepath = joinpath(testpath, "valvedata_wrappers_1.csv")
valves = read_valves(path = valvepath, delim = ',', skiplines = 1) #implicit read_valves test

model = WellModel(wellbore = well, roughness = 0.0, valves = valves,
                    temperature_method = "Shiu", geothermal_gradient = 1.0, BHT = 200,
                    pressurecorrelation = HagedornAndBrown, WHP = 350 - pressure_atmospheric, dp_est = 25,
                    q_o = 0, q_w = 500, GLR = 4500, APIoil = 35, sg_water = 1.0, sg_gas = 0.8, CHP = 1000, naturalGLR = 0)

tubing_pressures, casing_pressures, valvedata = gaslift_model!(model, find_injectionpoint = true, dp_min = 100) #also an implied test for 100% water cut

Δmds = [MDs[i] - MDs[i-1] for i in 81:length(MDs)]
ΔPs = [tubing_pressures[i] - tubing_pressures[i-1] for i in 81:length(MDs)]
gradients = ΔPs ./ Δmds
mean(x) = sum(x) / length(x)

expected_gradient = 0.433 / GouldWaterVolumeFactor(mean(tubing_pressures[81:end]), mean(model.temperatureprofile[81:end]))
@test mean(gradients) ≈ expected_gradient atol = 0.005

valve_table(valvedata)

#%% implicit plot test
if test_plots
    plot_gaslift(model.wellbore, tubing_pressures, casing_pressures, model.temperatureprofile, valvedata, nothing) |> x->draw(SVG("plot-gl-core.svg", 5inch, 4inch), x)
end

end #testset for gas lift wrappers
