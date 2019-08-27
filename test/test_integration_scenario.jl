# helper function for testing to load the example results
function load_example(path, ncol, delim = ',', skiplines = 1) #make sure to use nscenario+1 cols

    nlines = countlines(path) - skiplines
    output = Array{Float64, 2}(undef, nlines, ncol)
    filestream = open(path, "r")

    try
        for skip in 1:skiplines
            readline(filestream)
        end

        for (index, line) in enumerate(eachline(filestream))
            parsedline = parse.(Float64, split(line, delim, keepempty = false))
            output[index, :] = parsedline
        end
    finally
        close(filestream)
    end

    return output
end


# helper function for testing to interpolate the results to match wellbore segmentation
function match_example(well::Wellbore, example::Array{Float64,2})
    nrow = length(well.md)
    ncol = size(example)[2]
    output = Array{Float64, 2}(undef, nrow, ncol)
    output[:, 1] = well.md

    output[1, 2:end] = example[1, 2:end]
    for depth_index in 2:length(well.md)
        depth = well.md[depth_index]
        row_above = example[example[:,1] .<= depth, :][end,:]

        if row_above[1] == depth
            output[depth_index,2:end] = row_above[2:end]
        else
            row_below = example[example[:,1] .> depth, :][1,:]
            #interpolated_row = y1 + (y2 - y1)/(x2 - x1) * (depth - x1) :
            interpolated_values = [row_above[i] + (row_below[i] - row_above[i])/(row_below[1] - row_above[1]) * (depth - row_above[1]) for i in 2:ncol]
            output[depth_index,2:end] = interpolated_values
        end
    end

    return output
end


#%% general parameters
dp_est = 10 #psi
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


@testset "Beggs and Brill with Palmer correction - scenarios" begin

# load test results to compare & generate interpolations from expected results at same depths as test well segmentation
examplepath = joinpath(testpath,"Perform_BandB_with_Palmer_correction.csv")
test_example = load_example(examplepath, length(scenarios)+1)
matched_example = match_example(testwell, test_example)
matched_example[:,2:end] = matched_example[:,2:end] .- 14.65 #convert to psig

# generate test data -- map across all examples
corr = BeggsAndBrill
test_results = Array{Float64, 2}(undef, length(testwell.md), 1 + length(scenarios))
test_results[:,1] = testwell.md
for (index, scenario) in enumerate(scenarios)
    WC = parameters[:WC][scenario]
    rate = parameters[:rate][scenario]
    q_o = rate * (1 - WC)
    q_w = rate * WC
    GLR = parameters[:GLR][scenario]

    test_results[:,index+1] = traverse_topdown(wellbore = testwell, roughness = roughness, temperatureprofile = temp_profiles[scenario],
                             pressurecorrelation = corr, dp_est = dp_est, error_tolerance = error_tolerance,
                             q_o = q_o, q_w = q_w, GLR = GLR, APIoil = oil_API, sg_water = sg_water, sg_gas = sg_gas,
                             WHP = outlet_pressure, molFracCO2 = CO2)
end

# compare test data
# NOTE: points after the first survey point at > 90° inclination are discarded;
# the reference data appears to force positive frictional effects.
hz_index = findnext(x -> x >= 90, testwell.inc, 1)
compare_tolerance = 40 #every scenario except the high-rate scenario can take a 15-psi tolerance
for index in 2:length(scenarios)+1
    println("Comparing scenario ", index-1, " on index ", index)
    @test all( abs.(test_results[1:hz_index,index] .- matched_example[1:hz_index,index])  .<= compare_tolerance )
    println("Max difference : ", maximum(abs.(test_results[1:hz_index,index] .- matched_example[1:hz_index,index])))
end

end #testset for B&B scenarios


@testset "Hagedorn & Brown with G&W correction - scenarios" begin

# load test results to compare & generate interpolations from expected results at same depths as test well segmentation
examplepath = joinpath(testpath,"Perform_HandB_with_GriffithWallis.csv")
test_example = load_example(examplepath, length(scenarios)+1)
matched_example = match_example(testwell, test_example)
matched_example[:,2:end] = matched_example[:,2:end] .- 14.65 #convert to psig

# generate test data -- map across all examples
corr = HagedornAndBrown
test_results = Array{Float64, 2}(undef, length(testwell.md), 1 + length(scenarios))
test_results[:,1] = testwell.md
for (index, scenario) in enumerate(scenarios)
    WC = parameters[:WC][scenario]
    rate = parameters[:rate][scenario]
    q_o = rate * (1 - WC)
    q_w = rate * WC
    GLR = parameters[:GLR][scenario]

    test_results[:,index+1] = traverse_topdown(wellbore = testwell, roughness = roughness, temperatureprofile = temp_profiles[scenario],
                             pressurecorrelation = corr, dp_est = dp_est, error_tolerance = error_tolerance,
                             q_o = q_o, q_w = q_w, GLR = GLR, APIoil = oil_API, sg_water = sg_water, sg_gas = sg_gas,
                             WHP = outlet_pressure, molFracCO2 = CO2)
end

# compare test data
# NOTE: points after the first survey point at > 90° inclination are discarded
hz_index = findnext(x -> x >= 90, testwell.inc, 1)
compare_tolerance = 85 #lloser tolerance to H&B due to wide range of methods to calculate correlating groups & reynolds numbers
for index in 2:length(scenarios)+1
    println("Comparing scenario ", index-1, " on index ", index)
    @test all( abs.(test_results[1:hz_index,index] .- matched_example[1:hz_index,index])  .<= compare_tolerance )
    println("Max difference : ", maximum(abs.(test_results[1:hz_index,index] .- matched_example[1:hz_index,index])))
end

end #H&B testset
