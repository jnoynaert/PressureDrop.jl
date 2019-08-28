#test flags:

const devmode = false #test only subsets
const devtests = ("test_pvt.jl") #tuple of filenames to run for limited-subset tests
const test_plots = false

const run_benchmarks = false; timelimit = 5 #time limit in seconds for each benchmarking process

using Test
using PressureDrop

if test_plots #note that doc generation on deployment implicitly tests all plotting functions
    using Gadfly
end

if devmode
    include.(devtests);
else
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
end

if run_benchmarks
    using BenchmarkTools
    include("runbenchmarks.jl")
end