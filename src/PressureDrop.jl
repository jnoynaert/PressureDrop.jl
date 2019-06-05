# Package for computing multiphase pressure profiles for gas lift optimization of oil & gas wells

module PressureDrop

using Requires, PrettyTables

import Base.show #to export type printing methods

export  Wellbore, GasliftValves, WellModel,
        traverse_topdown, casing_traverse_topdown, read_survey, pressure_and_temp,
        plot_pressure, plot_pressures, plot_temperature, plot_pressureandtemp, plot_gaslift,
        valve_table, estimate_valve_Rvalue,
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
        ChenFrictionFactor


#%% setup
include("types.jl")
include("utilities.jl")
include("pvtproperties.jl")
include("valvecalculations.jl")
include("pressurecorrelations.jl")
include("tempcorrelations.jl")
include("casingcalculations.jl")


function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("plottingfunctions.jl")
end


#%% TODO: Don't kill traverses. Leave as is and make a "mapper wrapper" to convert struct to args for people who just want a traverse
#model struct. ONLY applies to wrapper functions.
"""
Makes it easier to iterate well models
"""
mutable struct WellModel

    wellbore::Wellbore; roughness
    valves::Union{GasliftValves, missing}
    temperatureprofile::Union{Array{Real,1}, missing}
    temperature_method; WHT; geothermal_gradient; BHT; casing_temp_factor
    pressurecorrelation::Function; outlet_referenced::Bool
    WHP; CHP; dp_est; dp_est_inj; error_tolerance; error_tolerance_inj
    q_o; q_w; GLR
    injection_point; naturalGLR
    APIoil; sg_water; sg_gas; sg_gas_inj
    molFracCO2; molFracH2S; molFracCO2_inj; molFracH2S_inj
    pseudocrit_pressure_correlation::Function
    pseudocrit_temp_correlation::Function
    Z_correlation::Function
    gas_viscosity_correlation::Function
    solutionGORcorrelation::Function
    oilVolumeFactor_correlation::Function
    waterVolumeFactor_correlation::Function
    dead_oil_viscosity_correlation::Function
    live_oil_viscosity_correlation::Function
    frictionfactor::Function

    function WellModel(;wellbore, roughness, valves = missing, temperatureprofile = missing,
                        temperature_method = "linear", WHT = missing, geothermal_gradient = missing, BHT = missing, casing_temp_factor = 0.85,
                        pressurecorrelation = BeggsAndBrill, outlet_referenced = true,
                        WHP, CHP = missing, dp_est, dp_est_inj = 0.1 * dp_est, error_tolerance = 0.1, error_tolerance_inj = 0.05,
                        q_o, q_w, GLR, injection_point = missing, naturalGLR = missing,
                        APIoil, sg_water, sg_gas, sg_gas_inj = sg_gas,
                        molFracCO2 = 0.0, molFracH2S = 0.0, molFracCO2_inj = molFracCO2, molFracH2S_inj = molFracH2S,
                        pseudocrit_pressure_correlation = HankinsonWithWichertPseudoCriticalPressure,
                        pseudocrit_temp_correlation = HankinsonWithWichertPseudoCriticalTemp,
                        Z_correlation = KareemEtAlZFactor, gas_viscosity_correlation = LeeGasViscosity,
                        solutionGORcorrelation = StandingSolutionGOR, oilVolumeFactor_correlation = StandingOilVolumeFactor,
                        waterVolumeFactor_correlation = GouldWaterVolumeFactor,
                        dead_oil_viscosity_correlation = GlasoDeadOilViscosity, live_oil_viscosity_correlation = ChewAndConnallySaturatedOilViscosity,
                        frictionfactor = SerghideFrictionFactor)

        new(wellbore, roughness, valves, temperatureprofile, temperature_method, WHT, geothermal_gradient, BHT,
            pressurecorrelation, outlet_referenced, casing_temp_factor, WHP, CHP, dp_est, dp_est_inj, error_tolerance, error_tolerance_inj,
            q_o, q_w, GLR, injection_point, naturalGLR, APIoil, sg_water, sg_gas, sg_gas_inj, molFracCO2, molFracH2S, molFracCO2_inj, molFracH2S_inj,
            pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation, gas_viscosity_correlation, solutionGORcorrelation,
            oilVolumeFactor_correlation, waterVolumeFactor_correlation, dead_oil_viscosity_correlation, live_oil_viscosity_correlation, frictionfactor)
    end

end

function Base.show(io::IO, model::WellModel)

    fields = fieldnames(WellModel)
    values = map(f -> getfield(model, f), fields)
    msg = String.(fields) .* " : " .* string.(values) .* "\n"

    print(io, msg)
end

macro @run(input, func)

    return :( begin
    fields = fieldnames(typeof($input))
    values = map(f -> getfield($input, f), fields)

    args = (;(f=>v for (f,v) in zip(fields,values))...) #splat and append to convert to NamedTuple

    $func(;args...)

    end)
end
#TODO: you need an @pressuredrop macro or similar that strips all the fields from a wellmodel struct and calls the function
#this will also prevent the requirement of rewriting all your existing tests? maybe doesn't need to be a macro.
#KEY poitn of this todo is to generate a kwargs-type namedtuple from the model and pass that to the traverse functions
#see https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
#except you can't do the above unless your function accepts extra kwargs...but maybe that is okay?
#yes, it is...and you can use that approach to every single named-arg to provide complete pass-through. brilliant. documentation remains the same.


#TODO: add testing for the above for a master wrapper


#%% core functions
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

    while abs(dp_est - dp_calc) > error_tolerance
        dp_est = dp_calc
        p_avg = p_initial + dp_est/2

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

Develop pressure traverse from wellhead down to datum in psia, returning a pressure profile as an Array{Float64,1}.

Pressure correlation functions available:
- `BeggsAndBrill` with Payne correction factors
- `HagedornAndBrown` with Griffith and Wallis bubble flow correction

# Arguments

All arguments are named keyword arguments.

## Required
- `wellbore::Wellbore`: Wellbore object that defines segmentation/mesh, with md, tvd, inclination, and hydraulic diameter
- `roughness`: pipe wall roughness in inches
- `temperatureprofile::Array{Float64, 1}`: temperature profile (in °F) as an array with **matching entries for each pipe segment defined in the Wellbore input**
- `outlet_pressure`: absolute outlet pressure (wellhead pressure) in **psia**
- `dp_est`: estimated starting pressure differential (in psi) to use for all segments--impacts convergence time
- `q_o`: oil rate in stocktank barrels/day
- `q_w`: water rate in stb/d
- `GLR`: **total** wellhead gas:liquid ratio, inclusive of injection gas, in scf/bbl
- `APIoil`: API gravity of the produced oil
- `sg_water`: specific gravity of produced water
- `sg_gas`: specific gravity of produced gas

## Optional
- `injection_point = missing`: injection point in MD for gas lift, above which total GLR is used, and below which natural GLR is used
- `naturalGLR = missing`: GLR to use below point of injection, in scf/bbl
- `pressurecorrelation::Function = BeggsAndBrill: pressure correlation to use
- `error_tolerance = 0.1`: error tolerance for each segment in psi
- `molFracCO2 = 0.0`, `molFracH2S = 0.0`: produced gas fractions of hydrogen sulfide and CO2, [0,1]
- `pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure`: psuedocritical pressure function to use
- `pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp`: pseudocritical temperature function to use
- `Z_correlation::Function = KareemEtAlZFactor`: natural gas compressibility/Z-factor correlation to use
- `gas_viscosity_correlation::Function = LeeGasViscosity`: gas viscosity correlation to use
- `solutionGORcorrelation::Function = StandingSolutionGOR`: solution GOR correlation to use
- `oilVolumeFactor_correlation::Function = StandingOilVolumeFactor`: oil volume factor correlation to use
- `waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor`: water volume factor correlation to use
- `dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity`: dead oil viscosity correlation to use
- `live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity`: saturated oil viscosity correction function to use
- `frictionfactor::Function = SerghideFrictionFactor`: correlation function for Darcy-Weisbach friction factor
"""
function traverse_topdown(;wellbore::Wellbore, roughness, temperatureprofile::Array{Float64, 1},
                            pressurecorrelation::Function = BeggsAndBrill,
                            outlet_pressure, dp_est, error_tolerance = 0.1,
                            q_o, q_w, GLR, injection_point = missing, naturalGLR = missing,
                            APIoil, sg_water, sg_gas, molFracCO2 = 0.0, molFracH2S = 0.0,
                            pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure, pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp,
                            Z_correlation::Function = KareemEtAlZFactor, gas_viscosity_correlation::Function = LeeGasViscosity, solutionGORcorrelation::Function = StandingSolutionGOR,
                            oilVolumeFactor_correlation::Function = StandingOilVolumeFactor, waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor,
                            dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity, live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity, frictionfactor::Function = SerghideFrictionFactor,
                            kwargs...) #catch extra arguments from a WellModel for convenience

    nsegments = length(wellbore.md)

    @assert nsegments == length(temperatureprofile) "Number of wellbore segments does not match number of temperature points."

    if !ismissing(injection_point) && !ismissing(naturalGLR)
        inj_index = searchsortedlast(wellbore.md, injection_point)

        if well.md[inj_index] != injection_point
            if (injection_point - wellbore.md[inj_index]) > (wellbore.md[inj_index+1] - injection_point) #choose closest point
                inj_index += 1
            end

            @info """Specified injection point at $injection_point not explicitly included in wellbore. Using $(well.md[inj_index]) as an approximate match.
            Use the Wellbore constructor with a set of gas lift valves to add precise injection points."""
        end

        GLRs = vcat(repeat([GLR], inner = inj_index), repeat([naturalGLR], inner = nsegments - inj_index))
    elseif !ismissing(injection_point) || !ismissing(naturalGLR)
        @info "Both an injection point and natural GLR should be specified--ignoring partial specification."
        GLRs = repeat([GLR], inner = nsegments)
    else #no injection point
        GLRs = repeat([GLR], inner = nsegments)
    end


    pressures = Array{Float64, 1}(undef, nsegments)
    pressure_initial = pressures[1] = outlet_pressure

    @inbounds for i in 2:nsegments
        dp_calc = calculate_pressuresegment_topdown(pressurecorrelation, pressure_initial, dp_est,
                                                    (temperatureprofile[i] + temperatureprofile[i-1])/2, #average temperature
                                                    wellbore.md[i-1], wellbore.md[i], wellbore.tvd[i-1], wellbore.tvd[i],
                                                    (wellbore.inc[i] + wellbore.inc[i-1])/2, #average inclination between survey points
                                                    wellbore.id[i], roughness,
                                                    q_o, q_w, GLRs[i], APIoil, sg_water, sg_gas, molFracCO2, molFracH2S,
                                                    pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation,
                                                    gas_viscosity_correlation, solutionGORcorrelation, oilVolumeFactor_correlation, waterVolumeFactor_correlation,
                                                    dead_oil_viscosity_correlation, live_oil_viscosity_correlation, frictionfactor, error_tolerance)

        pressure_initial += dp_calc
        pressures[i] = pressure_initial
    end

    return pressures
end


"""
traverse_topdown(;model::WellModel)

calculate top-down traverse from a WellModel object. Requires the following fields to be defined in the model:

...
"""
function traverse_topdown(model::WellModel)

    @run model traverse_topdown
end


"""
BHP_summary(pressures, well)

Print the summary for a bottomhole pressure traverse of a well.
"""
function BHP_summary(pressures, well)
    println("Flowing bottomhole pressure of $(round(pressures[end], digits = 1)) psia at $(well.md[end])' MD.",
        "\nAverage gradient $(round(pressures[end]/well.md[end], digits = 3)) psi/ft (MD), $(round(pressures[end]/well.tvd[end], digits = 3)) psi/ft (TVD).")
end


#TODO: update this documentation to reflect using a WellModel and mutating temps
"""
pressure_and_temp(;model::WellModel)

Develop pressure traverse in psia and temperature profile in °F from wellhead down to datum for a WellModel object. Requires the following fields to be defined in the model:

Returns a pressure profile as an Array{Float64,1} and a temperature profile as an Array{Float64,1}, referenced to the measured depths in the original Wellbore object.

Pressure correlation functions available:
- `BeggsAndBrill` with Payne correction factors
- `HagedornAndBrown` with Griffith and Wallis bubble flow correction

Temperature methods available:
- "Shiu" to utilize the Ramey 1962 method with the Shiu 1980 relaxation factor correlation
- "linear" for a linear interpolation between wellhead and bottomhole temperature based on TVD

# Arguments

All arguments are named keyword arguments.

## Required
- `well::Wellbore`: Wellbore object that defines segmentation/mesh, with md, tvd, inclination, and hydraulic diameter
- `roughness`: pipe wall roughness in inches
- `temperature_method = "linear"`: temperature method to use; "Shiu" for Ramey method with Shiu relaxation factor, "linear" for linear interpolation
- `WHT = missing`: wellhead temperature in °F; required for `temperature_method = "linear"`
- `geothermal_gradient = missing`: geothermal gradient in °F per 100 ft; required for `temperature_method = "Shiu"`
- `BHT` = bottomhole temperature in °F
- `WHP`: absolute outlet pressure (wellhead pressure) in **psia**
- `dp_est`: estimated starting pressure differential (in psi) to use for all segments--impacts convergence time
- `q_o`: oil rate in stocktank barrels/day
- `q_w`: water rate in stb/d
- `GLR`: **total** wellhead gas:liquid ratio, inclusive of injection gas, in scf/bbl
- `APIoil`: API gravity of the produced oil
- `sg_water`: specific gravity of produced water
- `sg_gas`: specific gravity of produced gas

## Optional
- `injection_point = missing`: injection point in MD for gas lift, above which total GLR is used, and below which natural GLR is used
- `naturalGLR = missing`: GLR to use below point of injection, in scf/bbl
- `pressurecorrelation::Function = BeggsAndBrill: pressure correlation to use
- `error_tolerance = 0.1`: error tolerance for each segment in psi
- `molFracCO2 = 0.0`, `molFracH2S = 0.0`: produced gas fractions of hydrogen sulfide and CO2, [0,1]
- `pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure`: psuedocritical pressure function to use
- `pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp`: pseudocritical temperature function to use
- `Z_correlation::Function = KareemEtAlZFactor`: natural gas compressibility/Z-factor correlation to use
- `gas_viscosity_correlation::Function = LeeGasViscosity`: gas viscosity correlation to use
- `solutionGORcorrelation::Function = StandingSolutionGOR`: solution GOR correlation to use
- `oilVolumeFactor_correlation::Function = StandingOilVolumeFactor`: oil volume factor correlation to use
- `waterVolumeFactor_correlation::Function = GouldWaterVolumeFactor`: water volume factor correlation to use
- `dead_oil_viscosity_correlation::Function = GlasoDeadOilViscosity`: dead oil viscosity correlation to use
- `live_oil_viscosity_correlation::Function = ChewAndConnallySaturatedOilViscosity`: saturated oil viscosity correction function to use
- `frictionfactor::Function = SerghideFrictionFactor`: correlation function for Darcy-Weisbach friction factor
- `outlet_referenced = true`: whether to use outlet pressure (WHP) or inlet pressure (BHP) for
"""
function pressure_and_temp!(m::WellModel)

    if m.temperature_method == "linear"
        @assert !(any(ismissing.((m.WHT,m.BHT)))) "Must specific a wellhead temperature & BHT to utilize linear temperature method."
        m.temperatureprofile = @run m linear_wellboretemp
    elseif m.temperature_method == "Shiu"
        @assert !(any(ismissing.((BHT,geothermal_gradient)))) "Must specify a geothermal gradient & BHT to utilize Shiu/Ramey temperature method.\nRefer to published geothermal gradient maps for your region to establish a sensible default."
        m.temperatureprofile = @run m Shiu_wellboretemp
    else
        throw(ArgumentError("Invalid temperature method. Use one of (\"Shiu\", \"linear\")."))
    end

    pressures = @run m traverse_topdown
    BHP_summary(pressures, well)

    return pressures
end
#TODO: update test


"""
"""
function pressures_and_temp(m::WellModel)

    tubing_pressures = @run m traverse_topdown
    casing_pressures = casing_traverse_topdown(m)

    return tubing_pressures, casing_pressures
end
#TODO: add test


#TODO: docs
"""
"""
function gaslift_model!(m::WellModel, find_injectionpoint::Bool = false, dp_min = 100)

    tubing_pressures = @run m pressure_and_temp!
    casing_pressures = casing_traverse_topdown(m)
    valvedata, injection_depth = valve_calcs(m.valves, m.well, m.sg_gas_inj, tubing_pressures, casing_pressures, m.temperatureprofile, m.temperatureprofile .* casing_temp_factor,
                            dp_min)

    if find_injectionpoint
        m.injection_point = injection_depth
        tubing_pressures = @run m traverse_topdown
        casing_pressures = casing_traverse_topdown(m)
        valvedata, _ = valve_calcs(m.valves, m.well, m.sg_gas_inj, tubing_pressures, casing_pressures, m.temperatureprofile, m.temperatureprofile .* casing_temp_factor,
                                dp_min)
    end

    return tubing_pressures, casing_pressures, valvedata
end
#TODO: test


end #module PressureDrop
