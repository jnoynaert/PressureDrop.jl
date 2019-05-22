"""
minimum implementation

Assumes no friction loss

Assumes no liquid entrained in gas stream
"""
function casing_pressure_segment_topdown(p_initial, dp_est, t_avg,
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

    while abs(dp_est - dp_calc) >= error_tolerance
        dp_est = dp_calc
        p_avg = p_initial + dp_est/2

        Z = Z_correlation(P_pc, T_pc, p_avg, t_avg)
        ρ_g = gasDensity_insitu(sg_gas, Z, p_avg, t_avg)

        dp_calc = (1/144.0) * ρ_g * dh_tvd
    end

    return dp_calc
end

# add friction version via multiple dispatch by adding necessary parameters
