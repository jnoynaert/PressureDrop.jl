using Test

include("pvptproperties.jl")

#%% Gas

@test LeeGasViscosity(0.7, 2000, 160, 0.75) ≈ 0.017 atol = 0.001

@test HankinsonWithWichertPseudoCriticalTemp(0.65, 0.05, 0.08) .- (370.2, 351.4, 18.09) |> x -> all(abs.(x) .<= 1)

@test HankinsonWithWichertPseudoCriticalPressure(0.65, 0.05, 0.08) ≈ 634 atol = 2

@test PapayZFactor

@test gasVolumeFactor

@test gasDensity

#%% Oil
@test StandingSolutionGOR

@test StandingOilVolumeFactor

@test oilDensity

@test BeggsAndRobinsonDeadOilViscosity

@test GlasoDeadOilViscosity

@test ChewAndConnallySaturatedOilViscosity

#%% Water

@test waterDensity_stb

@test GouldWaterVolumeFactor

@test assumedWaterViscosity == 1.0

@test waterDensity
