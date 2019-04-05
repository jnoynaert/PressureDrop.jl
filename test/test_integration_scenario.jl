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


#=General parameters
BHT = 165

WHP = 220 PSIG

oil_API = 35
sg_gas = 0.65
CO2 = 0.5%
N2 = 2%
sg_water = 1.07
roughness = 0.0006500

settings for PVT in IHS:
dead oil visco -> Glaso
Saturated oil vis -> chew & ChewAndConnallySaturatedOilViscosity
undersaturated -> vazquez
gas -> Lee
Water -> matthews & russell
oil density -> standing
bpp & R_S -> standing
oil comp -> vazquez & beggs
oil fvf -> standing
Z -> Hall & Yarborough
=

survey = Sawgrass 9 to TD (12175 MD), first segment edited to 0 inclination (460 md/tvd)

tubing id = 2.441
=#
#%% general parameters
BHT = 165

#%% scenarios
scenarios = (:A, :B, :C, :D, :E)

parameters =    (rate = (A = 500, B = 250, C = 1000, D = 3000, E = 50),
                 WC = (A = 0.5, B = 0.25, C = 0.75, D = 0.85, E = 0.25),
                 GLR = (A = 4500, B = 6000, C = 3000, D = 1200, E = 10000),
                 WHT = (A = 100, B = 90, C = 105, D = 115, E = 80))

#%% load test Wellbore
testpath = "../test/testdata/Sawgrass 9-32/"
surveypath = testpath*"Test survey - sawgrass 9.csv"
testwell = read_survey(path = surveypath, id = 2.441)

#%% generate temperature profiles
function create_temps(scenario)
        WHT = parameters[:WHT][scenario]
        temp_slope = (BHT - WHT) / maximum(testwell.tvd)
        return [WHT + depth * temp_slope for depth in testwell.tvd]
end

temps = [create_temps(s) for s in scenarios]
temp_profiles = NamedTuple{scenarios}(temps)


#%% example: test scenario A for B&B
#TODO: finisht test:
# load test results to compare & generate interpolations from expected results at same depths as test well segmentation
examplepath = testpath*"Perform - B&B with Palmer correction.csv"
BB_example = load_example(examplepath, length(scenarios)+1)
BB_matched_example = match_example(testwell, BB_example)

# generate test data -- map across all examples
#traverse_topdown(wellbore = testwell, roughness = 0.0006, temperatureprofile = test_temp,
                             pressurecorrelation = BeggsAndBrill, dp_est = 10, error_tolerance = 0.1,
                             q_o = 400, q_w = 500, GLR = 2000, APIoil = 36, sg_water = 1.05, sg_gas = 0.75,
                             outlet_pressure = 150)

# compare test data
tolerance = 35 #psi
#abs.(results .- example) .<= tolerance
