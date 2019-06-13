"""
calculate_casing_pressuresegment_topdown(<arguments>)

Helper function to calculate the pressure drop for a single casing segment containing only injection gas, using a inlet-referenced (casing head referenced) approach.

Assumes no friction or entrained liquid -- uses density only.

Pressure inputs are in **psia**.

See `casing_traverse_topdown`.
"""
function casing_pressuresegment_topdown(p_initial, dp_est, t_avg,
                                            tvd_initial, tvd_end,
                                            sg_gas, molFracCO2, molFracH2S,
                                            pseudocrit_pressure_correlation::Function, pseudocrit_temp_correlation::Function, Z_correlation::Function,
                                            error_tolerance = 0.1)

    dh_tvd = tvd_end - tvd_initial
    p_avg = p_initial + dp_est/2

    P_pc = pseudocrit_pressure_correlation(sg_gas, molFracCO2, molFracH2S)
    _, T_pc, _ = pseudocrit_temp_correlation(sg_gas, molFracCO2, molFracH2S)
    Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
    ρ_g = gasDensity_insitu(sg_gas, Z, p_avg, t_avg)

    dp_calc = (1/144.0) * ρ_g * dh_tvd

    while abs(dp_est - dp_calc) > error_tolerance
        dp_est = dp_calc
        p_avg = p_initial + dp_est/2

        Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
        ρ_g = gasDensity_insitu(sg_gas, Z, p_avg, t_avg)

        dp_calc = (1/144.0) * ρ_g * dh_tvd
    end

    return dp_calc
end


"""
casing_traverse_topdown(;<named arguments>)

Develops pressure traverse from casing head down to datum in psia, returning a pressure profile as an Array{Float64,1}.

Uses only density and is only applicable to pure gas injection, i.e. assumes no friction loss and no liquid entrained in gas stream (reasonable assumptions for relatively dry gas taken through several compression stages and injected through relatively large casing).

Pressure inputs are in **psig**.

# Arguments

All arguments are named keyword arguments.

## Required
- `wellbore::Wellbore`: Wellbore object that defines segmentation/mesh, with md, tvd, inclination, and hydraulic diameter
- `temperatureprofile::Array{Float64, 1}`: temperature profile (in °F) as an array with **matching entries for each pipe segment defined in the Wellbore input**
- `CHP`: casing head pressure, i.e. absolute surface injection pressure in **psig**
- `dp_est`: estimated starting pressure differential (in psi) to use for all segments--impacts convergence time
- `sg_gas`: specific gravity of produced gas

## Optional
- `error_tolerance = 0.1`: error tolerance for each segment in psi
- `molFracCO2 = 0.0`, `molFracH2S = 0.0`: produced gas fractions of hydrogen sulfide and CO2, [0,1]
- `pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure`: psuedocritical pressure function to use
- `pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp`: pseudocritical temperature function to use
- `Z_correlation::Function = KareemEtAlZFactor`: natural gas compressibility/Z-factor correlation to use
"""
function casing_traverse_topdown(;wellbore::Wellbore, temperatureprofile::Array{Float64, 1},
                                    CHP, dp_est, error_tolerance = 0.1,
                                    sg_gas, molFracCO2 = 0.0, molFracH2S = 0.0,
                                    pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure, pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp,
                                    Z_correlation::Function = KareemEtAlZFactor)

    CHP += pressure_atmospheric

    nsegments = length(wellbore.md)

    @assert nsegments == length(temperatureprofile) "Number of wellbore segments does not match number of temperature points."

    pressures = Array{Float64, 1}(undef, nsegments)
    pressure_initial = pressures[1] = CHP

    @inbounds for i in 2:nsegments
        dp_calc = casing_pressuresegment_topdown(pressure_initial, dp_est,
                                                    (temperatureprofile[i] + temperatureprofile[i-1])/2, #average temperature
                                                    wellbore.tvd[i-1], wellbore.tvd[i],
                                                    sg_gas, molFracCO2, molFracH2S,
                                                    pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation,
                                                    error_tolerance)

        pressure_initial += dp_calc
        pressures[i] = pressure_initial
    end

    return pressures .- pressure_atmospheric
end


"""
Remaps casing traverse to work with WellModels
"""
function casing_traverse_topdown(m::WellModel)

    casing_traverse_topdown(;wellbore = m.wellbore, temperatureprofile = m.temperatureprofile .* m.casing_temp_factor,
                                        CHP = m.CHP, dp_est = m.dp_est_inj, error_tolerance = m.error_tolerance_inj,
                                        sg_gas = m.sg_gas_inj, molFracCO2 = m.molFracCO2_inj, molFracH2S = m.molFracH2S_inj,
                                        pseudocrit_pressure_correlation = m.pseudocrit_pressure_correlation,
                                        pseudocrit_temp_correlation = m.pseudocrit_temp_correlation,
                                        Z_correlation = m.Z_correlation)
end
