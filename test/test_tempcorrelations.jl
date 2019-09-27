include("../src/pressurecorrelations.jl")
include("../src/tempcorrelations.jl")

@testset "Temperature correlations" begin

q_o = 2219
q_w = 11
APIoil = 24
sg_water = 1
GLR = 1.762 * 1000^2 / (2219 + 11)
sg_gas = 0.7
id = 2.992
whp = 150 #psig
A = Shiu_Beggs_relaxationfactor(q_o, q_w, GLR, APIoil, sg_water, sg_gas, id, whp)

@test A ≈ 2089 atol = 1

g_g = 0.0106 * 100
bht = 173
z = 5985 - 0

@test Ramey_temp(z, bht, A, g_g) ≈ 130 atol = 2
@test Ramey_temp(1, bht, A, g_g) ≈ 173 atol = 0.1

#= visual test
depths = range(1, stop = z, length = 100) ;

test_temps = Ramey_wellboretemp.(depths, 0, bht, A, g_g)
depths_plot = z .- depths |> collect

using Gadfly

plot(y = depths_plot, x = test_temps, Geom.line, Coord.cartesian(yflip = true))
=#

end #Temperature correlations
