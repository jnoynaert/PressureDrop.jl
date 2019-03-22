using Test
using PressureDrop

include("test_pvt.jl")
include("test_pressurecorrelations.jl")
include("test_tempcorrelations.jl")
include("test_full.jl") #tolerance test
#include("test_integration.jl")

#TODO: change all of your test function calls for unexported functions to be fully qualifeid,
# and remove the file imports.
