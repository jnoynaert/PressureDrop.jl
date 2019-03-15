module PressureDrop

# Calculate pressure and temperature profiles for oil & gas wells.

#TODO: avoid duplicating calculations
#TODO: Finish profiling using @time, @code_warntype, Traceur, ProfileView
#TODO: @inbounds master loop
#TODO: add test/runtests.jl per Takacs results, or per a reference result from IHS, and run on every new build once you have a main loop in place

#using
#import


include("pvtproperties.jl")
include("pressurecorrelations.jl")
include("tempcorrelations.jl")



struct wellbore
    md::Array{Float64, 1}
    tvd::Array{Float64, 1}
    inc::Array{Float64, 1}
    id::Array{Float64, 1}

    function wellbore(md, tvd, inc, id::Array{Float64, 1})
        lens = length.([md, tvd, inc, id])

        return count(x -> x == lens[1], lens) == length(lens) ?
            new(md, tvd, inc, id) :
            throw(DimensionMismatch("Mismatched number of wellbore elements used in wellbore constructor."))
    end
end

wellbore(md, tvd, inc, id::Float64) = wellbore(md, tvd, inc, repeat([id], inner = length(md))) #convenience constructor for uniform tubulars

test() = return 1, 2, 3


"""
helper fn for calculating pressure segment

assumes you are using helper fns of the top-down form that return ΔP and not dpdh
"""
function calculate_pressuresegment_topdown(pressurecorrelation, p_initial, dp_est, t_avg,
                                            md_initial, md_end, tvd_initial, tvd_end, id, args...)

    dh_md = md_end - md_initial
    dh_tvd = tvd_end - tvd_initial
    p_avg = p_initial + dp_est/2
    uphill_flow = inclination <= 90.0 ? true : false

    #calculate pvt properties at p_avg, t_avg FROM PVT FNS PASSED AS ARGS :


    P_pc = HankinsonWithWichertPseudoCriticalPressure(sg_gas, molFracCO2, molFracH2S)
    _, T_pc, _ = HankinsonWithWichertPseudoCriticalTemp(sg_gas, molFracCO2, molFracH2S)
    Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
    B_g = gasVolumeFactor(p_avg, Z, t_avg)
    μ_g = gas_viscosity_correlation(sg_gas, p_avg, t_avg, Z)
    R_s = solutionGORcorrelation(...)
    v_sg = gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)
    B_o = oilVolumeFactor_correlation(args...)
    B_w = waterVolumeFactor_correlation(args...)
    v_sl = liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)
    ρ_l = mixture_properties_simple(q_o, q_w, oilDensity_insitu(...), waterDensity_insitu(...))
    σ_l = mixture_properties_simple(q_o, q_w, gas_oil_interfacialtension(APIoil, p_avg, t_avg), gas_water_interfacialtension(pressureAbsolute, tempF))
    μ_l = mixture_properties_simple(q_o, q_w, oil_viscosity_correlation(...), assumedWaterViscosity) #may need to use a composed function rather than separate live/dead oil versions

    dp_calc = pressurecorrelation(dh_md, dh_tvd, inclination, id,
                                    v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g, roughness, p_avg,
                                    uphill_flow)

    while abs(dp_est - dp_calc) >= error_tolerance:
        dp_est = dp_calc
        p_avg = p_initial + dp_est/2

        #TODO: recalc PRESSURE DEPENDENT PROPERTIES ONLY...except that looks like...all of them?

        dp_calc = pressurecorrelation(args..., dp_est)
    end

    return max(dp_calc, 0) #constrain top-down drop for producers in the initial case
end



"""
Develop pressure traverse from wellhead down to datum.

This function does not handle re-meshing/segmenting, calculating BHT, etc
    ^^will outsource and create a higher level wrapper if needed

Inputs: wellbore: Array{Float64, segments, 4}: [md_array, tvd_array, inc_array, id_array],
outlet_pressure, temp profile,
produced GOR,
gas s.g., water s.g., oil API
which correlations to use for fluid properties
which pressure drop correlation to use as the function name

Returns: pressure profile, temperature profile as two separate Array{Float64, 1}s.
"""
function traverse_topdown(;survey, pressurecorrelation = BeggsAndBrill, dp_est = bestguessforpressurestep)
    # TODO: use GOR to initialize model?, since GLR is ambiguous term for gas lift opt, but GOR is understand to be reservoir #s
    # TODO: make sure to utilize named arguments

    nsegments = length(wellbore.md)

    @assert nsegments == length(temperatureprofile) "Number of wellbore segments does not match number of temperature segments."

    pressures = Array{Float64, 1}(undef, nsegments)
    pressure_initial = pressures[1] = outlet_pressure

    for i in 2:nsegments
        #TODO: be very, very sure the arguments match here--you want to use positional for performance but could easily get a nigh-untrackable error
        #TODO: go back and validate where you are returning absolute pressures vs ΔPs, and where you are picking avg or point pressure/temp spots
        dp_calc = calculate_pressuresegment_topdown(pressurecorrelation, pressure_initial, dp_est, temperatureprofile[i],
                                                            wellbore.md[i-1], wellbore.md[i], wellbore.tvd[i-1], wellbore.tvd[i], (wellbore.id[i] + wellbore.id[i-1])/2,
                                                            args...)
        pressure_initial += dp_calc
        pressures[i] = pressure_initial
    end

    return pressures
end

end #module PressureDrop
