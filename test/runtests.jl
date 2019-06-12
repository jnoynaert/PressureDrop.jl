using Test
using PressureDrop

const test_plots = false

if test_plots
    using Gadfly
end

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
