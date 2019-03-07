# Pressure correlations for PressureDrop package.

#TODO: ensure appropriate kwarg handling so that all correlations utilize a common interface.


"""
Beggs and Brill flow pattern as a string.

Takes no-slip holdup (λ_l) and mixture Froude number (N_Fr).

Beggs and Brill. Takacs p87.
"""
function BeggsAndBrillFlowMap(λ_l, N_Fr) #graphical test bypassed in test suite--rerun if modifying this function

    if N_Fr < 316 * λ_l ^ 0.302 && N_Fr < 9.25e-4 * λ_l^-2.468
        return "segrated"
    elseif N_Fr >= 9.25e-4 * λ_l^-2.468 && N_Fr < 0.1 * λ_l^-1.452
        return "transition"
    elseif N_Fr >= 316 * λ_l ^ 0.302 || N_Fr >= 0.5 * λ_l^-6.738
        return "distributed"
    else
        return "intermittent"
    end
end


"""
Documentation here.

Note you will use modified version from Takacs p89.
"""
function BeggsAndBrill()

    g = 32.2 #gravitational constant ft/s
    α = (1-inclination) * π / 180 #inclination in rad measured from horizontal

    λ_l = v_sl / v_m #no-slip liquid holdup
    N_Fr = v_m^2/(g * d) #mixture Froude number #d is pipe diameter in inches

    flowpattern = BeggsAndBrillFlowMap(λ_l, N_Fr)

    dpdl_el = (1/144.0) * ρ_m * sin(α) #elevation component
    dpdl_f = 1.294e-3 * f * (ρ_ns * v_m^2) / d #frictional component
    E_k = 2.16e-4 * f * (v_m * v_sg * ρ_ns) / p #kinetic energy/acceleration component

    dp_dl = (dpdl_el*d + dpdl_f*d) / (1-E_k)

    return dp_dl
end
