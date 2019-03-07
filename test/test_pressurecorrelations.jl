using Test

include("pressurecorrelations.jl")

#%% Beggs and Brill flowmap
@test BeggsAndBrillFlowMap(x, y) == "segregated"
@test BeggsAndBrillFlowMap(x, y) == "transition"
@test BeggsAndBrillFlowMap(x, y) == "intermittent"
@test BeggsAndBrillFlowMap(x, y) == "distributed"

#= B&B flowmap graphical test, to make sure the correct mapping is generated: passed 3/7/2019.
# Re-run after any modifications to B&B flowmap function.
# Note that this test does not check for stability at the region boundaries.

λ_l = 10 .^ collect(-4:0.1:0);
N_Fr = 10 .^ collect(-1:0.1:3);

function crossjoin(x, y)

    output = Array{Float64, 2}(undef, length(x) * length(y), 2)

    index = 1
    for i in x
        for j in y
            output[index, 1] = i
            output[index, 2] = j
            index += 1
        end
    end

    return output
end

testdata = crossjoin(λ_l, N_Fr);
flowpattern = map(BeggsAndBrillFlowMap, testdata[:,1], testdata[:,2]);


using Gadfly

plot(x = testdata[:,1], y = testdata[:,2], color = flowpattern,
Geom.rectbin,
Scale.x_log10, Scale.y_log10, Coord.Cartesian(xmin=-4, xmax=0, ymin=-1, ymax=3))
=#

#%% Beggs and Brill simplified case

# test case conditions for Takacs un-inclined example (example parameters from pg 47)
inclination = 90.0
md = 8245.0
tvd = 8245.0
q_oil = 375.0
GOR = 480.0
oil_sg = 0.82
gas_sg = 0.916
WHP = 346.6
BHP = 2529.5
WHT = 80.3
BHT = 140.3

# Standing correlations at surface
oil_API = 41
R_s = 109.7
B_o = 1.05
A_pipe = 0.0325

v_sl = 0.79

p_pc = 655.8
t_pc = 451.9

p_pr = 0.528
t_pr = 1.19

Z = 0.883
B_g = 0.0389
v_sg = 1.93

#Standing correlations at bottomhole
R_s = 480
B_o = 1.286
v_sl = 0.96

@test
