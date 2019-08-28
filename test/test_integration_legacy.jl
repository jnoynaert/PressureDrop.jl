@testset "Takacs B&B vertical single large segment" begin

pressurecorrelation = BeggsAndBrill
p_initial = 346.6 - 14.65
dp_est = 100
t_avg = 107.2
md_initial = 0
md_end = 3700
tvd_initial = 0
tvd_end = 3700
inclination = 0
id = 2.441
roughness = 0.0003
q_o = 375
q_w = 0
GLR = 480
APIoil = 41.06
sg_water = 1
sg_gas = 0.916
molFracCO2 = 0
molFracH2S = 0
pseudocrit_pressure_correlation = HankinsonWithWichertPseudoCriticalPressure
pseudocrit_temp_correlation = HankinsonWithWichertPseudoCriticalTemp
Z_correlation = PapayZFactor
gas_viscosity_correlation = LeeGasViscosity
solutionGORcorrelation = StandingSolutionGOR
bubblepoint = 3000
oilVolumeFactor_correlation = StandingOilVolumeFactor
waterVolumeFactor_correlation = GouldWaterVolumeFactor
dead_oil_viscosity_correlation = GlasoDeadOilViscosity
live_oil_viscosity_correlation = ChewAndConnallySaturatedOilViscosity
frictionfactor = SerghideFrictionFactor
error_tolerance = 1.0 #psi

P_pc = pseudocrit_pressure_correlation(sg_gas, molFracCO2, molFracH2S)
    _, T_pc, _ = pseudocrit_temp_correlation(sg_gas, molFracCO2, molFracH2S)

ΔP_est = PressureDrop.calculate_pressuresegment(pressurecorrelation, p_initial, dp_est, t_avg,
                                md_end - md_initial, tvd_end - tvd_initial, inclination, true, id, roughness,
                                q_o, q_w, GLR, GLR, APIoil, sg_water, sg_gas,
                                Z_correlation, P_pc, T_pc,
                                gas_viscosity_correlation, solutionGORcorrelation, bubblepoint, oilVolumeFactor_correlation, waterVolumeFactor_correlation,
                                dead_oil_viscosity_correlation, live_oil_viscosity_correlation, frictionfactor, error_tolerance)

@test ΔP_est ≈ 3770 * (0.170 + 0.312) / 2 atol = 100 #order of magnitude test based on average between example points in Takacs (50).

end #testset Takacs B&B vertical



@testset "IHS Cleveland 6 - B&B full wellbore" begin

#%% end to end test

testpath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata/Cleveland_6/Test_survey_Cleveland_6.csv")

testwell = read_survey(path = testpath, id_included = false, maxdepth = 10000, id = 2.441)
test_temp = collect(range(85, 160, length = length(testwell.md)))

pressure_values = traverse_topdown(wellbore = testwell, roughness = 0.0006, temperatureprofile = test_temp,
                                    pressurecorrelation = BeggsAndBrill, dp_est = 10, error_tolerance = 0.1,
                                    q_o = 400, q_w = 500, GLR = 2000, APIoil = 36, sg_water = 1.05, sg_gas = 0.75,
                                    WHP = 150 - 14.65)

ihs_data = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata/Cleveland_6/Perform_results_cleveland_6_long.csv")
ihs_pressures = [parse.(Float64, split(line, ',', keepempty = false)) for line in readlines(ihs_data)[2:end]] |>
                        x -> hcat(x...)'


@test (ihs_pressures[end, 2] - 14.65) ≈ pressure_values[end] atol = 25

model = WellModel(wellbore = testwell, roughness = 0.0006, temperatureprofile = test_temp,
                    pressurecorrelation = BeggsAndBrill, dp_est = 10, error_tolerance = 0.1,
                    q_o = 400, q_w = 500, GLR = 2000, APIoil = 36, sg_water = 1.05, sg_gas = 0.75,
                    WHP = 150 - 14.65)

@test (ihs_pressures[end, 2] - 14.65) ≈ traverse_topdown(model)[end] atol = 25

#= view results
ihs_temps = "C:/pressuredrop.git/test/testdata/Cleveland_6/Perform_temps_cleveland_6_long.csv"
ihs_temps = [parse.(Float64, split(line, ',', keepempty = false)) for line in readlines(ihs_temps)[2:end]] |>
                        x -> hcat(x...)' ;

using Gadfly
set_default_plot_size(8.5inch, 11inch)

plot(   layer(x = pressure_values, y = testwell.md, Geom.line),
        layer(x = test_temp, y = testwell.md, Geom.line, Theme(default_color = "red")),
        layer(x = ihs_pressures[:,2], y = ihs_pressures[:,1], Geom.line, Theme(default_color = "green")),
        layer(x = ihs_temps[:,2], y = ihs_temps[:,1], Geom.line, Theme(default_color = "orange")),
        Coord.cartesian(yflip = true)) #matches


[pressure_values[i] - pressure_values[i-1] for i in 2:length(pressure_values)] #negative drops only occur where inclination is negative
=#

end #testset IHS Cleveland 6 - B&B
