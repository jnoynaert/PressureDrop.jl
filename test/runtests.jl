#test flags:

devmode = false #test only subsets
isCI = get(ENV, "CI", nothing) == "true" #check if running in Travis
devtests = ("test_pvt.jl") #tuple of filenames to run for limited-subset tests
test_plots = false #falls through to test_wrappers.jl

run_benchmarks = false; timelimit = 5 #time limit in seconds for each benchmarking process

using Test
using PressureDrop

if test_plots #note that doc generation on deployment implicitly tests all plotting functions
    using Gadfly
end

if isCI || !devmode
    include("test_types.jl")
    include("test_utilities.jl")
    include("test_pvt.jl")
    include("test_pressurecorrelations.jl")
    include("test_tempcorrelations.jl")
    include("test_casingcalcs.jl")
    include("test_valvecalcs.jl")
    include("test_integration_legacy.jl") #tolerance test
    include("test_integration_scenario.jl")
    include("test_wrappers.jl")
    include("test_regressions.jl")
 else #devmode
    include.(devtests);
end

if run_benchmarks
    using BenchmarkTools
    include("runbenchmarks.jl")
end
