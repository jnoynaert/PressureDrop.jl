# run benchmarks on integration scenarios

BenchmarkTools.DEFAULT_PARAMETERS.seconds = timelimit
Base.CoreLogging.disable_logging(Base.CoreLogging.Info) #remove info-level outputs to avoid console spam
@warn "Setting up integration benchmarks..."

#%% general parameters
dp_est = 10. #psi
error_tolerance = 0.1 #psi
outlet_pressure = 220 - 14.65 #WHP in psig
oil_API = 35
sg_gas = 0.65
CO2 = 0.005
sg_water = 1.07
roughness = 0.0006500
BHT = 165
id = 2.441

#%% scenario parameters
scenarios = (:A, :B, :C, :D, :E)

parameters =    (rate = (A = 500, B = 250, C = 1000, D = 3000, E = 50),
                 WC = (A = 0.5, B = 0.25, C = 0.75, D = 0.85, E = 0.25),
                 GLR = (A = 4500, B = 6000, C = 3000, D = 1200, E = 10000),
                 WHT = (A = 100, B = 90, C = 105, D = 115, E = 80))

#%% load test Wellbore
testpath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata/Sawgrass_9_32")
surveypath = joinpath(testpath, "Test_survey_Sawgrass_9.csv")
testwell = read_survey(path = surveypath, id = id)

#%% generate linear temperature profiles
function create_temps(scenario)
        WHT = parameters[:WHT][scenario]

        return linear_wellboretemp(WHT = WHT, BHT = BHT, wellbore = testwell)
end

temps = [create_temps(s) for s in scenarios]
temp_profiles = NamedTuple{scenarios}(temps) #temp profiles labelled by scenario

@warn "Benchmarking Beggs & Brill..."
corr = BeggsAndBrill
timings = Array{Float64,1}(undef, length(scenarios))
for (index, scenario) in enumerate(scenarios)
    WC = parameters[:WC][scenario]
    rate = parameters[:rate][scenario]
    q_o = rate * (1 - WC)
    q_w = rate * WC
    GLR = parameters[:GLR][scenario]
    temps = temp_profiles[scenario]

    timings[index] = @belapsed begin traverse_topdown(wellbore = testwell, roughness = roughness, temperatureprofile = $temps,
                                pressurecorrelation = corr, dp_est = dp_est, error_tolerance = error_tolerance,
                                q_o = $q_o, q_w = $q_w, GLR = $GLR, APIoil = oil_API, sg_water = sg_water, sg_gas = sg_gas,
                                WHP = outlet_pressure, molFracCO2 = CO2)
                        end

end

@warn "$corr timing | min: $(minimum(timings)) s | max: $(maximum(timings)) s | mean: $(sum(timings)/length(timings)) s"


@warn "Benchmarking Hagedorn & Brown..."
corr = HagedornAndBrown
timings = Array{Float64,1}(undef, length(scenarios))
for (index, scenario) in enumerate(scenarios)
    WC = parameters[:WC][scenario]
    rate = parameters[:rate][scenario]
    q_o = rate * (1 - WC)
    q_w = rate * WC
    GLR = parameters[:GLR][scenario]
    temps = temp_profiles[scenario]

    timings[index] = @belapsed begin traverse_topdown(wellbore = testwell, roughness = roughness, temperatureprofile = $temps,
                                pressurecorrelation = corr, dp_est = dp_est, error_tolerance = error_tolerance,
                                q_o = $q_o, q_w = $q_w, GLR = $GLR, APIoil = oil_API, sg_water = sg_water, sg_gas = sg_gas,
                                WHP = outlet_pressure, molFracCO2 = CO2)
                        end
end

@warn "$corr timing | min: $(minimum(timings)) s | max: $(maximum(timings)) s | mean: $(sum(timings)/length(timings)) s"