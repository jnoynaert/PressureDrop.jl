# Package for computing multiphase pressure profiles for gas lift optimization of oil & gas wells

module PressureDrop

push!(LOAD_PATH, @__DIR__) #enable separate loading of PressurePlots.jl

import Base.show

export  Wellbore, traverse_topdown, read_survey,
        BeggsAndBrill,
        HagedornAndBrown,
        Shiu_wellboretemp, Ramey_temp, Shiu_Beggs_relaxationfactor, linear_wellboretemp,
        LeeGasViscosity,
        HankinsonWithWichertPseudoCriticalTemp,
        HankinsonWithWichertPseudoCriticalPressure,
        PapayZFactor,
        KareemEtAlZFactor,
        KareemEtAlZFactor_simplified,
        StandingSolutionGOR,
        StandingOilVolumeFactor,
        BeggsAndRobinsonDeadOilViscosity,
        GlasoDeadOilViscosity,
        ChewAndConnallySaturatedOilViscosity,
        GouldWaterVolumeFactor,
        SerghideFrictionFactor,
        ChenFrictionFactor#,
        #TODO: export other wrappers


"""
struct Wellbore: object to define a flow path as an input for pressure drop calculations

See `read_survey` for helper method to create a Wellbore object from deviation survey files.

# Fields
-`md::Array{Float64, 1}`: measured depth for each segment in feet
-`inc::Array{Float64, 1}`: inclination from vertical for each segment in degrees, e.g. true vertical = 0°
-`tvd::Array{Float64, 1}`: true vertical depth for each segment in feet
-`id::Array{Float64, 1}`: inner diameter for each pip segment in inches

# Constructors
`Wellbore(md, inc, tvd, id::Array{Float64, 1})`: defines a new Wellbore object from a survey with inner diameter defined for each segment. Lengths of each input array must be equal.

`Wellbore(md, inc, tvd, id::Float64)`: defines a new Wellbore object with a uniform ID along the entire flow path.
"""
struct Wellbore
    md::Array{Float64, 1}
    inc::Array{Float64, 1}
    tvd::Array{Float64, 1}
    id::Array{Float64, 1}

    function Wellbore(md, inc, tvd, id::Array{Float64, 1})
        lens = length.([md, inc, tvd, id])

        return count(x -> x == lens[1], lens) == length(lens) ?
            new(md, inc, tvd, id) :
            throw(DimensionMismatch("Mismatched number of wellbore elements used in wellbore constructor."))
    end
end


#convenience constructor for uniform tubulars
Wellbore(md, inc, tvd, id::Float64) = Wellbore(md, inc, tvd, repeat([id], inner = length(md)))

#Printing for Wellbore structs
Base.show(io::IO, well::Wellbore) = print(io,
    "Wellbore with $(length(well.md)) segments.\nEnds at $(well.md[end])' MD / $(well.tvd[end])' TVD. \nMax inclination $(maximum(well.inc))°. Average ID $(round(sum(well.id)/length(well.id), digits = 3)) in.")

include("pvtproperties.jl")
include("pressurecorrelations.jl")
include("tempcorrelations.jl")
include("utilities.jl")


"""
calculate_pressuresegment_topdown(<arguments>)

Helper function to calculate the pressure drop for a single pipe segment, using an outlet-referenced approach.

Method:
1. Generates PVT properties for the average conditions in the segment for a given estimate for the pressure drop in the second.
2. Calculates a new estimated pressure drop using the PVT properties.
3. Compares the original and new pressure drop estimates to validate the stability of the estimate.
4. Iterate through steps 1-3 until a stable estimate is found (indicated by obtaining a difference between pre- and post-PVT that is within the given error tolerance).

See traverse_topdown for a full enumeration of arguments.
"""
function calculate_pressuresegment_topdown(pressurecorrelation::Function, p_initial, dp_est, t_avg,
                                            md_initial, md_end, tvd_initial, tvd_end, inclination, id, roughness,
                                            q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2, molFracH2S,
                                            pseudocrit_pressure_correlation::Function, pseudocrit_temp_correlation::Function, Z_correlation::Function,
                                            gas_viscosity_correlation::Function, solutionGORcorrelation::Function, oilVolumeFactor_correlation::Function, waterVolumeFactor_correlation::Function,
                                            dead_oil_viscosity_correlation::Function, live_oil_viscosity_correlation::Function, frictionfactor::Function, error_tolerance = 0.1)

    dh_md = md_end - md_initial
    dh_tvd = tvd_end - tvd_initial
    p_avg = p_initial + dp_est/2
    uphill_flow = inclination <= 90.0

    P_pc = pseudocrit_pressure_correlation(sg_gas, molFracCO2, molFracH2S)
    _, T_pc, _ = pseudocrit_temp_correlation(sg_gas, molFracCO2, molFracH2S)
    Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
    ρ_g = gasDensity_insitu(sg_gas, Z, p_avg, t_avg)
    B_g = gasVolumeFactor(p_avg, Z, t_avg)
    μ_g = gas_viscosity_correlation(sg_gas, p_avg, t_avg, Z)
    R_s = solutionGORcorrelation(APIoil, sg_gas, p_avg, t_avg)
    v_sg = gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)
    B_o = oilVolumeFactor_correlation(APIoil, sg_gas, R_s, p_avg, t_avg)
    B_w = waterVolumeFactor_correlation(p_avg, t_avg)
    v_sl = liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)
    ρ_l = mixture_properties_simple(q_o, q_w, oilDensity_insitu(APIoil,  sg_gas,  R_s,  B_o), waterDensity_insitu(sg_water, B_w))
    σ_l = mixture_properties_simple(q_o, q_w, gas_oil_interfacialtension(APIoil, p_avg, t_avg), gas_water_interfacialtension(p_avg, t_avg))
    μ_oD = dead_oil_viscosity_correlation(APIoil, t_avg)
    μ_l = mixture_properties_simple(q_o, q_w, live_oil_viscosity_correlation(μ_oD, R_s), assumedWaterViscosity)

    dp_calc = pressurecorrelation(dh_md, dh_tvd, inclination, id,
                                    v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g, roughness, p_avg, frictionfactor,
                                    uphill_flow)

    while abs(dp_est - dp_calc) >= error_tolerance
        dp_est = dp_calc
        p_avg = p_initial + dp_est/2

        P_pc = pseudocrit_pressure_correlation(sg_gas, molFracCO2, molFracH2S)
        _, T_pc, _ = pseudocrit_temp_correlation(sg_gas, molFracCO2, molFracH2S)
        Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
        ρ_g = gasDensity_insitu(sg_gas, Z, p_avg, t_avg)
        B_g = gasVolumeFactor(p_avg, Z, t_avg)
        μ_g = gas_viscosity_correlation(sg_gas, p_avg, t_avg, Z)
        R_s = solutionGORcorrelation(APIoil, sg_gas, p_avg, t_avg)
        v_sg = gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)
        B_o = oilVolumeFactor_correlation(APIoil, sg_gas, R_s, p_avg, t_avg)
        B_w = waterVolumeFactor_correlation(p_avg, t_avg)
        v_sl = liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)
        ρ_l = mixture_properties_simple(q_o, q_w, oilDensity_insitu(APIoil,  sg_gas,  R_s,  B_o), waterDensity_insitu(sg_water, B_w))
        σ_l = mixture_properties_simple(q_o, q_w, gas_oil_interfacialtension(APIoil, p_avg, t_avg), gas_water_interfacialtension(p_avg, t_avg))
        μ_oD = dead_oil_viscosity_correlation(APIoil, t_avg)
        μ_l = mixture_properties_simple(q_o, q_w, live_oil_viscosity_correlation(μ_oD, R_s), assumedWaterViscosity)

        dp_calc = pressurecorrelation(dh_md, dh_tvd, inclination, id,
                                        v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g, roughness, p_avg, frictionfactor,
                                        uphill_flow)
    end

    return dp_calc #allows negatives
end


"""
traverse_topdown(;<named arguments>)

Develop pressure traverse from wellhead down to datum, returning a pressure profile as an Array{Float64,1}.

# Arguments

All arguments are named keyword arguments.

## Required
-`wellbore::Wellbore`: Wellbore object that defines segmentation/mesh, with md, tvd, inclination, and hydraulic diameter
-`roughness`: pipe wall roughness in inches
-`temperatureprofile::Array{Float64, 1}`: temperature profile (in °F) as an array with **matching entries for each pipe segment defined in the Wellbore input**
-`outlet_pressure`: outlet pressure (wellhead pressure) in psia
-`dp_est`: estimated starting pressure differential (in psi) to use for all segments--impacts convergence time
-`q_o`: oil rate in stocktank barrels/day
-`q_w`: water rate in stb/d
-`GLR`: producing gas:liquid ratio at standard conditions in scf/bbl
-`APIoil`: API gravity of the produced oil
-`sg_water`: specific gravity of produced water
-`sg_gas`: specific gravity of produced gas

## Optional
-`pressurecorrelation::Function: pressure correlation to use
-`error_tolerance = 0.1`: error tolerance for each segment in psi
-`molFracCO2 = 0.0`, `molFracH2S = 0.0`: produced gas fractions of hydrogen sulfide and CO2, [0,1]
-`pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure`: psuedocritical pressure function to use
-`pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp`: pseudocritical temperature function to use
-`Z_correlation::Function = KareemEtAlZFactor`: natural gas compressibility/Z-factor correlation to use
-`gas_viscosity_correlation::Function = LeeGasViscosity`: gas viscosity correlation to use
-`solutionGORcorrelation::Function = StandingSolutionGOR`: solution GOR correlation to use
-`oilVolumeFactor_correlation::Function = StandingOilVolumeFactor`: oil volume factor correlation to use
-`waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor`: water volume factor correlation to use
-`dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity`: dead oil viscosity correlation to use
-`live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity`: saturated oil viscosity correction function to use
-`frictionfactor::Function = SerghideFrictionFactor`: correlation function for Darcy-Weisbach friction factor
"""
function traverse_topdown(;wellbore::Wellbore, roughness, temperatureprofile::Array{Float64, 1},
                            pressurecorrelation::Function = BeggsAndBrill,
                            outlet_pressure, dp_est, error_tolerance = 0.1,
                            q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2 = 0.0, molFracH2S = 0.0,
                            pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure, pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp,
                            Z_correlation::Function = KareemEtAlZFactor, gas_viscosity_correlation::Function = LeeGasViscosity, solutionGORcorrelation::Function = StandingSolutionGOR,
                            oilVolumeFactor_correlation::Function = StandingOilVolumeFactor, waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor,
                            dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity, live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity, frictionfactor::Function = SerghideFrictionFactor)

    nsegments = length(wellbore.md)

    @assert nsegments == length(temperatureprofile) "Number of wellbore segments does not match number of temperature points."

    pressures = Array{Float64, 1}(undef, nsegments)
    pressure_initial = pressures[1] = outlet_pressure

    @inbounds for i in 2:nsegments
        dp_calc = calculate_pressuresegment_topdown(pressurecorrelation, pressure_initial, dp_est,
                                                    (temperatureprofile[i] + temperatureprofile[i-1])/2, #average temperature
                                                    wellbore.md[i-1], wellbore.md[i], wellbore.tvd[i-1], wellbore.tvd[i],
                                                    (wellbore.inc[i] + wellbore.inc[i-1])/2, #average inclination between survey points
                                                    wellbore.id[i], roughness,
                                                    q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2, molFracH2S,
                                                    pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation,
                                                    gas_viscosity_correlation, solutionGORcorrelation, oilVolumeFactor_correlation, waterVolumeFactor_correlation,
                                                    dead_oil_viscosity_correlation, live_oil_viscosity_correlation, frictionfactor, error_tolerance)

        pressure_initial += dp_calc
        pressures[i] = pressure_initial
    end

    return pressures
end


#TODO: pressure + temp wrapper function from Wellbore with passthrough arguments; allow linear or Shiu method temp


end #module PressureDrop
