"""
Assumes flow has stabilized and transient time component does not change.
WHP in psiA
GLR in scf/bbl

surprisingly sensitive to gas rate
"""
function Shiu_Beggs_relaxationfactor(q_o, q_w, APIoil, sg_water, GLR, sg_gas, id, whp)
    outlet_pressure = whp - 14.67
    sg_oil = 141.5/(APIoil + 131.5)
    q_g = (q_o + q_w) * GLR
    w = 1/86400 * (350*(q_o*sg_oil + q_w*sg_water) + 0.0764*q_g*sg_gas) #mass flow rate in lb/sec
    ρ_l_sc = 62.4 * mixture_properties_simple(q_o, q_w, sg_oil, sg_water)

    # Original Shiu-Beggs coefficients for known WHP:
    C_0 = 0.0063
    C_1 = 0.4882
    C_2 = 2.9150
    C_3 = -0.3476
    C_4 = 0.2219
    C_5 = 0.2519
    C_6 = 4.7240

    return C_0 * w^C_1 * ρ_l_sc^C_2 * id^C_3 * outlet_pressure^C_4 * APIoil^C_5 * sg_gas^C_6 #relaxation distance
end


"""
geothermal gradient, g_g (°F per 100 ft)
tvd FROM well bottom, z (ft)
relaxation factor, A
"""
function Ramey_wellboretemp(z, inclination, T_bh, A, G_g = 1.0)

    g_g = G_g / 100
    α = (90 - inclination) * π / 180

    return T_bh - g_g * z * sin(α) + A * g_g * sin(α) * (1 - exp(-z / A))
end


"""
"""
function linear_wellboretemp(;WHT, BHT, well::Wellbore)
    temp_slope = (BHT - WHT) / maximum(well.tvd)

    return [WHT + depth * temp_slope for depth in well.tvd]
end


#TODO: function to create temp traverse in main file, and use a Wellbore struct as an input
"""
"""
