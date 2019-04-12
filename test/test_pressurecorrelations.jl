include("../src/pressurecorrelations.jl")

@testset "Friction factors" begin
    @test ChenFrictionFactor(35700, 2.259, 0.001) ≈ 0.0063 atol = 0.001
end

@testset "Superficial velocities" begin

    @test liquidvelocity_superficial(375, 0, 2.441, 1.05, 0) ≈ 0.79 atol = 0.01
    @test gasvelocity_superficial(375, 0, 480, 109.7, 2.441, 0.0389) ≈ 1.93 atol = 0.01

end

#%% Beggs and Brill flowmap

#=
@test BeggsAndBrillFlowMap(x, y) == "segregated"
@test BeggsAndBrillFlowMap(x, y) == "transition"
@test BeggsAndBrillFlowMap(x, y) == "intermittent"
@test BeggsAndBrillFlowMap(x, y) == "distributed"
=#

#= B&B flowmap graphical test, to make sure the correct mapping is generated: passed 3/7/2019.
# Re-run after any modifications to B&B flowmap function.
# Note that this does not establish stability at the region boundaries.

λ_l = 10 .^ collect(-4:0.1:0);
N_Fr = 10 .^ collect(-1:0.1:3);

function crossjoin(x, y)

    output = Array{Float64, 2}(undef, length(x) * length(y), 2)

    index = 1
    @inbounds for i in x
        @inbounds for j in y
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

@testset "Beggs and Brill vertical - single segment" begin

    inclination = 0
    v_sl = 0.79
    v_sg = 1.93
    v_m = 2.72
    ρ_l = 49.9
    ρ_g = 1.79
    id = 2.441
    σ_l = 8
    μ_l = 3
    μ_g = 0.02
    roughness = 0.0006
    pressure_est = 346.6

    @test BeggsAndBrill(1, 1, inclination, id, v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g, roughness, pressure_est, true, false) ≈ 0.170 atol = 0.01

end
