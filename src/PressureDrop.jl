module PressureDrop

push!(LOAD_PATH, @__DIR__) #enable separate loading of PressurePlots.jl

#TODO: avoid duplicating calculations
#TODO: add bubblepoint logic w/ appropriate multiphase calculation & undersaturated oil properties
#TODO: Finish profiling using @time, @code_warntype, Traceur, ProfileView
#TODO: @inbounds master loop
#TODO: add test/runtests.jl per Takacs results, or per a reference result from IHS, and run on every new build once you have a main loop in place

export  Wellbore, traverse_topdown, read_survey,
        BeggsAndBrill,
        HagedornAndBrown,
        Ramey_wellboretemp, Shiu_Beggs_relaxationfactor, linear_wellboretemp,
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
        ChenFrictionFactor


#TODO: pretty printing for Wellbore
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

Wellbore(md, inc, tvd, id::Float64) = Wellbore(md, inc, tvd, repeat([id], inner = length(md))) #convenience constructor for uniform tubulars


include("pvtproperties.jl")
include("pressurecorrelations.jl")
include("tempcorrelations.jl")
include("utilities.jl")


"""
helper fn for calculating pressure segment

assumes you are using helper fns of the top-down form that return ΔP and not dpdh
"""
function calculate_pressuresegment_topdown(pressurecorrelation::Function, p_initial, dp_est, t_avg,
                                            md_initial, md_end, tvd_initial, tvd_end, inclination, id, roughness,
                                            q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2, molFracH2S,
                                            pseudocrit_pressure_correlation::Function, pseudocrit_temp_correlation::Function, Z_correlation::Function,
                                            gas_viscosity_correlation::Function, solutionGORcorrelation::Function, oilVolumeFactor_correlation::Function, waterVolumeFactor_correlation::Function,
                                            dead_oil_viscosity_correlation::Function, live_oil_viscosity_correlation::Function, frictionfactor::Function, error_tolerance = 0.01)

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

    return dp_calc #allows negatives...
end



"""
Develop pressure traverse from wellhead down to datum.

This function does not handle re-meshing/segmenting, calculating BHT, etc
    ^^will outsource and create a higher level wrapper if needed

Inputs: wellbore: wellbore struct,
outlet_pressure, temp profile,
produced GOR,
gas s.g., water s.g., oil API
which correlations to use for fluid properties
which pressure drop correlation to use as the function name

Returns: pressure profile, temperature profile as two separate Array{Float64, 1}s.
"""
function traverse_topdown(;wellbore::Wellbore, roughness, temperatureprofile::Array{Float64, 1}, pressurecorrelation::Function = BeggsAndBrill, outlet_pressure, dp_est, error_tolerance = 0.1,
                            q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2 = 0.0, molFracH2S = 0.0,
                            pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure, pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp,
                            Z_correlation::Function = KareemEtAlZFactor, gas_viscosity_correlation::Function = LeeGasViscosity, solutionGORcorrelation::Function = StandingSolutionGOR,
                            oilVolumeFactor_correlation::Function = StandingOilVolumeFactor, waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor,
                            dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity, live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity, frictionfactor::Function = SerghideFrictionFactor)

    nsegments = length(wellbore.md)

    @assert nsegments == length(temperatureprofile) "Number of wellbore segments does not match number of temperature segments."

    pressures = Array{Float64, 1}(undef, nsegments)
    pressure_initial = pressures[1] = outlet_pressure

    @inbounds for i in 2:nsegments
        dp_calc = calculate_pressuresegment_topdown(pressurecorrelation, pressure_initial, dp_est, temperatureprofile[i],
                                                    wellbore.md[i-1], wellbore.md[i], wellbore.tvd[i-1], wellbore.tvd[i], (wellbore.inc[i] + wellbore.inc[i-1])/2,
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


end #module PressureDrop
