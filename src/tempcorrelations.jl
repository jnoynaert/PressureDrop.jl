# temperature correlation functions for PressureDrop package.

"""
`Shiu_Beggs_relaxationfactor(<arguments>)`

Generates the relaxation factor, A, needed for the Ramey method, for underspecified conditions.

This correlation assumes flow has stabilized and that the transient time component f(t) is not changing.

# Arguments
All arguments are in U.S. field units.

- `q_o`: oil rate in stb/d
- `q_w`: water rate in stb/d
- `APIoil`: API oil gravity
- `sg_water`: water specific gravity
- `GLR`: gas:liquid ratio in scf/stb
- `sg_gas`: gas specific gravity
- `id`: flow path inner diameter in inches
- `WHP`: wellhead/outlet absolute pressure in **psig**
"""
function Shiu_Beggs_relaxationfactor(q_o, q_w, GLR, APIoil, sg_water, sg_gas, id, WHP)

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

    return C_0 * w^C_1 * ρ_l_sc^C_2 * id^C_3 * WHP^C_4 * APIoil^C_5 * sg_gas^C_6 #relaxation distance
end


"""
`Ramey_wellboretemp(z, inclination, T_bh, A, G_g = 1.0)`

Estimates wellbore temp using Ramey 1962 method.

# Arguments
- `z`: true vertical depth **from the bottom of the well**, ft
- `T_bh`: bottomhole temperature, °F
- `A`: relaxation factor
- `G_g = 1.0`: geothermal gradient in °F per 100 ft of true vertical depth
"""
function Ramey_temp(z, T_bh, A, G_g = 1.0)

    g_g = G_g / 100

    return T_bh - g_g * z + A * g_g * (1 - exp(-z / A))
end


"""
`linear_wellboretemp(;WHT, BHT, wellbore::Wellbore)`

Linear temperature profile from a wellhead temperature and bottomhole temperature in °F for a Wellbore object.

Interpolation is based on true vertical depth of the wellbore, not md.
"""
function linear_wellboretemp(;WHT, BHT, wellbore::Wellbore,
                            kwargs...) #catch extra arguments from a WellModel for convenience

    temp_slope = (BHT - WHT) / maximum(wellbore.tvd)

    return [WHT + depth * temp_slope for depth in wellbore.tvd]
end


"""
`Shiu_wellboretemp(<named arguments>)`

Wrapper to compute temperature profile for a Wellbore object using Ramey correlation with Shiu relaxation factor correlation.

# Arguments
- `BHT`: bottomhole temperature in °F
- `geothermal_gradient = 1.0`: geothermal gradient in °F per 100 feet
- `wellbore::Wellbore`: Wellbore object to use as reference for segmentation, inclination, and
- `q_o`: oil rate in stb/d
- `q_w`: water rate in stb/d
- `GLR`: gas:liquid ratio in scf/day
- `APIoil`: oil gravity
- `sg_water`: water specific gravity
- `sg_gas`: gas specific gravity
- `WHP`: wellhead/outlet absolute pressure in **psig**
"""
function Shiu_wellboretemp(;BHT, geothermal_gradient = 1.0, wellbore::Wellbore, q_o, q_w, GLR, APIoil, sg_water, sg_gas, WHP,
                            kwargs...) #catch extra arguments from a WellModel for convenience

    id_avg = sum(wellbore.id)/length(wellbore.id)
    A = Shiu_Beggs_relaxationfactor(q_o, q_w, GLR, APIoil, sg_water, sg_gas, id_avg, WHP) #use average inner diameter to calculate relaxation factor
    TD = maximum(wellbore.tvd)
    depths = TD .- wellbore.tvd
    temp_profile = [Ramey_temp(z, BHT, A, geothermal_gradient) for z in depths]

    return temp_profile
end
