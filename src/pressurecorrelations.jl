# Pressure correlations for PressureDrop package.


#%% Helper functions

"""
liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)

Returns superficial liquid velocity, v_sg.

Takes oil rate (q_o, stb/d), water rate (q_w, stb/d), pipe inner diameter (inches), oil volume factor (B_o, dimensionless) and water volume factor (B_w, dimensionless).

Note that this does not account for slip between liquid phases.
"""
function liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)
    A = π * (id/24.0)^2 #convert id in inches to ft

    if q_o > 0
        WOR = q_w / q_o
        return 6.5e-5 * (q_o + q_w) / A * (B_o/(1 + WOR) + B_w * WOR / (1 + WOR))
    else #100% WC
        return 6.5e-5 * q_w * B_w / A
    end
end


"""
gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)

Returns superficial liquid velocity, v_sg.

Takes oil rate, (q_o, stb/d), water rate (q_w, stb/d), gas:liquid ratio (scf/stb), solution gas:oil ratio (scf/stb), pipe inner diameter (inches), and gas volume factor (B_g).

Note that this does not account for slip between liquid phases.

"""
function gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)
    A = π * (id/24.0)^2 #convert id in inches to ft

    if q_o > 0
        WOR = q_w / q_o
        return 1.16e-5 * (q_o + q_w) / A * (GLR - R_s /(1 + WOR)) * B_g
    else #100% WC
        return 1.16e-5 * q_w * (GLR - R_s) * B_g / A
    end
end

# mixture velocity: just v_sg + v_sl


"""
mixture_properties_simple(q_o, q_w, property_o, property_w)

Weighted average for mixture properties.

Takes the oil and water rates in equivalent units, as well as their relative properties in equivalent units.

Does not account for oil slip, mixing effects, fluid expansion, behavior of emulsions, etc.
"""
function mixture_properties_simple(q_o, q_w, property_o, property_w)

    return (q_o * property_o + q_w * property_w) / (q_o + q_w)
end


"""
ChenFrictionFactor(N_Re, id, roughness = 0.01)

Uses the direct Chen 1979 correlation to determine friction factor, in place of the Colebrook implicit solution.

Takes the dimensionless Reynolds number, pipe inner diameter in *inches*, and roughness in *inches*.

Not intended for Reynolds numbers between 2000-4000.
"""
function ChenFrictionFactor(N_Re, id, roughness = 0.01) #Takacs #returns Economides value * 4

    if N_Re <= 4000 #laminar flow boundary ~2000-2300

        return 64 / N_Re
    else #turbulent flow
        k = roughness/id

        x = -2 * log10(k/3.7065 - 5.0452/N_Re * log10(k^1.1098 / 2.8257 + (7.149/N_Re)^0.8981))

        return 1/x^2
    end
end


"""
SerghideFrictionFactor(N_Re, id, roughness = 0.01)

Uses the direct Serghide 1984 correlation to determine friction factor, in place of the Colebrook implicit solution.

Takes the dimensionless Reynolds number, pipe inner diameter in *inches*, and roughness in *inches*.

Not intended for Reynolds numbers between 2000-4000.
"""
function SerghideFrictionFactor(N_Re, id, roughness = 0.01)
    if N_Re <= 4000 #laminar flow boundary ~2000-2300

        return 64 / N_Re
    else #turbulent flow
        k = roughness/id
        A = -2 * log10(k/3.7 + 12/N_Re)
        B = -2 * log10(k/3.7 + 2.51*A/N_Re)
        C = -2 * log10(k/3.7 + 2.51*B/N_Re)

        return (A-((B-A)^2 / (C - (2*B) + A)))^-2
    end
end




#%% Beggs and Brill

"""
BeggsAndBrillFlowMap(λ_l, N_Fr)

Beggs and Brill flow pattern as a string ∈ {"segregated", "transition", "distributed", "intermittent"}.

Takes no-slip holdup (λ_l) and mixture Froude number (N_Fr).
"""
function BeggsAndBrillFlowMap(λ_l, N_Fr) #graphical test bypassed in test suite--rerun if modifying this function

    if N_Fr < 316 * λ_l ^ 0.302 && N_Fr < 9.25e-4 * λ_l^-2.468
        return "segregated"
    elseif N_Fr >= 9.25e-4 * λ_l^-2.468 && N_Fr < 0.1 * λ_l^-1.452
        return "transition"
    elseif N_Fr >= 316 * λ_l ^ 0.302 || N_Fr >= 0.5 * λ_l^-6.738
        return "distributed"
    else
        return "intermittent"
    end
end


const BB_coefficients =     (segregated = (a = 0.980, b= 0.4846, c = 0.0868, e = 0.011, f = -3.7680, g = 3.5390, h = -1.6140),
                            intermittent = (a = 0.845, b = 0.5351, c = 0.0173, e = 2.960, f = 0.3050, g = -0.4473, h = 0.0978),
                            distributed = (a = 1.065, b = 0.5824, c = 0.0609),
                            downhill = (e = 4.700, f = -0.3692, g = 0.1244, h = -0.5056) )


"""
BeggsAndBrillAdjustedLiquidHoldup(flowpattern, λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection = true)

Helper function for Beggs and Brill. Returns adjusted liquid holdup, ε_l(α), with optional Payne et al correction applied to output.

Takes flow pattern (string ∈ {"segregated", "intermittent", "distributed"}), no-slip holdup (λ_l), Froude number (N_Fr),
liquid velocity number (N_lv), angle from horizontal (α, radians), uphill flow (boolean).
"""
function BeggsAndBrillAdjustedLiquidHoldup(flowpattern, λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection = true)

    if PayneCorrection && uphill_flow
        correctionfactor = 0.924
    elseif PayneCorrection
        correctionfactor = 0.685
    else
        correctionfactor = 1.0
    end

    flow = Symbol(flowpattern)
    a = BB_coefficients[flow][:a]
    b = BB_coefficients[flow][:b]
    c = BB_coefficients[flow][:c]

    ε_l_horizontal = a * λ_l^b / N_Fr^c #liquid holdup assuming horizontal (α = 0 rad)
    ε_l_horizontal = max(ε_l_horizontal, λ_l)

    if α ≈ 0 #horizontal flow
        return ε_l_horizontal
    else #inclined or vertical flow
        if uphill_flow
            if flowpattern == "distributed"
                ψ = 1.0
            else
                e = BB_coefficients[flow][:e]
                f = BB_coefficients[flow][:f]
                g = BB_coefficients[flow][:g]
                h = BB_coefficients[flow][:h]

                C = max( (1 - λ_l) * log(e * λ_l^f * N_lv^g * N_Fr^h), 0)

                if inclination ≈ 0 #vertical flow
                    ψ = 1 + 0.3 * C
                else
                    ψ = 1 + C * (sin(1.8*α) - (1/3) * sin(1.8*α)^3)
                end
            end
        else #downhill flow
            e = BB_coefficients[:downhill][:e]
            f = BB_coefficients[:downhill][:f]
            g = BB_coefficients[:downhill][:g]
            h = BB_coefficients[:downhill][:h]

            C = max( (1 - λ_l) * log(e * λ_l^f * N_lv^g * N_Fr^h), 0)

            if inclination ≈ 0 #vertical flow
                ψ = 1 + 0.3 * C
            else
                ψ = 1 + C * (sin(1.8*α) - (1/3) * sin(1.8*α)^3)
            end
        end

        return ε_l_horizontal * ψ * correctionfactor
    end

end



"""
BeggsAndBrill(<arguments>)

Calculates pressure drop for a single pipe segment using Beggs and Brill 1973 method, with optional Payne corrections.

Returns a ΔP in psi.

Doesn't account for oil/water phase slip, but does properly account for inclination.

As of release v0.9, assumes **outlet-defined** models only, i.e. top-down from wellhead; thus, uphill flow corresponds to producers and downhill flow to injectors.

For more information, see *Petroleum Production Systems* by Economides et al., or the the Fekete [reference on pressure drops](http://www.fekete.com/san/webhelp/feketeharmony/harmony_webhelp/content/html_files/reference_material/Calculations_and_Correlations/Pressure_Loss_Calculations.htm).


# Arguments
All arguments take U.S. field units.

- `md`: measured depth of the pipe segment, feet
- `tvd`: true vertical depth, feet
- `inclination`: inclination from vertical, degrees (e.g. vertical => 0)
- `id`: inner diameter of the pipe segment, inches
- `v_sl`: superficial liquid mixture velocity, ft/s
- `v_sg`: superficial gas velocity, ft/s
- `ρ_l`: liquid mixture density, lb/ft³
- `ρ_g`: gas density, lb/ft³
- `σ_l`: liquid/gas interfacial tension, centipoise
- `μ_l`: liquid mixture dynamic viscosity
- `μ_g`: gas dynamic viscosity
- `roughness`: pipe roughness, inches
- `pressure_est`: estimated average pressure of the pipe segment (needed to determine the kinetic effects component of the pressure drop)
- `frictionfactor::Function = SerghideFrictionFactor`: function used to determine the Moody friction factor
- `uphill_flow = true`: indicates uphill or downhill flow. It is assumed that the start of the 1D segment is an outlet and not an inlet
- `PayneCorrection = true`: indicates whether the Payne et al. 1979 corrections should be applied to prevent overprediction of liquid holdup.
"""
function BeggsAndBrill( md, tvd, inclination, id,
                        v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g,
                        roughness, pressure_est, frictionfactor::Function = SerghideFrictionFactor,
                        uphill_flow = true, PayneCorrection = true)

    α = (90 - inclination) * π / 180 #inclination in rad measured from horizontal

    #%% flow pattern and holdup:
    v_m = v_sl + v_sg
    λ_l = v_sl / v_m #no-slip liquid holdup
    N_Fr = 0.373 * v_m^2 / id #mixture Froude number #id is pipe diameter in inches
    N_lv = 1.938 * v_sl * (ρ_l / σ_l)^0.25 #liquid velocity number per Duns & Ros

    flowpattern = BeggsAndBrillFlowMap(λ_l, N_Fr)

    if flowpattern == "transition"
        B = (0.1 * λ_l^-1.4516 - N_Fr) / (0.1 * λ_l^-1.4516 - 9.25e-4 * λ_l^-2.468)
        ε_l_seg = BeggsAndBrillAdjustedLiquidHoldup("segregated", λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection)
        ε_l_int = BeggsAndBrillAdjustedLiquidHoldup("intermittent", λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection)
        ε_l_adj = B * ε_l_seg + (1 - B) * ε_l_int
    else
        ε_l_adj = BeggsAndBrillAdjustedLiquidHoldup(flowpattern, λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection)
    end

    if uphill_flow
        ε_l_adj = max(ε_l_adj, λ_l) #correction to original: for uphill flow, true holdup must by definition be >= no-slip holdup
    end # note that Payne factors reduce the overpredicted liquid holdups from the uncorrected form

    #%% friction factor:
    y = λ_l / ε_l_adj^2
    if 1.0 < y < 1.2
        s = log(2.2y - 1.2) #handle the discontinuity
    else
        ln_y = log(y)
        s = ln_y / (-0.0523 + 3.182 * ln_y - 0.872 * ln_y^2 + 0.01853 * ln_y^4)
    end

    fbyfn = exp(s) #normalizing friction factor f/fₙ

    ρ_ns = ρ_l * λ_l + ρ_g * (1-λ_l) #no-slip density
    μ_ns = μ_l * λ_l + μ_g * (1-λ_l) #no-slip friction in centipoise
    N_Re = 124 * ρ_ns * v_m * id / μ_ns #Reynolds number
    f_n = frictionfactor(N_Re, id, roughness)
    fric = f_n * fbyfn #friction factor


    #%% core calculation:
    ρ_m = ρ_l * ε_l_adj + ρ_g * (1 - ε_l_adj) #mixture density in lb/ft³

    dpdl_el = (1/144.0) * ρ_m
    friction_effect = uphill_flow ? 1 : -1 #note that friction MUST act against the direction of flow
    dpdl_f = friction_effect * 1.294e-3 * fric * (ρ_ns * v_m^2) / id #frictional component
    E_k = 2.16e-4 * fric * (v_m * v_sg * ρ_ns) / pressure_est #kinetic effects; typically negligible

    dp_dl = (dpdl_el * tvd + dpdl_f * md) / (1 - friction_effect*E_k) #assumes friction and kinetic effects both increase pressure in the same 1D direction

    return dp_dl
end #Beggs and Brill



#%% Hagedorn & Brown
"""
HagedornAndBrownLiquidHoldup(pressure_est, id, v_sl, v_sg, ρ_l, μ_l, σ_l)

Helper function to determine liquid holdup using the Hagedorn & Brown method.

Does not account for inclination or oil/water slip.
"""
function HagedornAndBrownLiquidHoldup(pressure_est, id, v_sl, v_sg, ρ_l, μ_l, σ_l)
    N_lv = 1.938 * v_sl * (ρ_l / σ_l)^0.25 #liquid velocity number per Duns & Ros
    N_gv = 1.938 * v_sg * (ρ_l / σ_l)^0.25 #gas velocity number per Duns & Ros; yes, use liquid density & viscosity
    N_d = 120.872 * id/12 * (ρ_l / σ_l)^0.5 #pipe diameter number; uses id in ft
    N_l = 0.15726 * μ_l * (1/(ρ_l * σ_l^3))^0.25 #liquid viscosity number

    CN_l = 0.061 * N_l^3 - 0.0929 * N_l^2 + 0.0505 * N_l + 0.0019 #liquid viscosity coefficient * liquid viscosity number

    H = N_lv / N_gv^0.575 * (pressure_est/14.7)^0.1 * CN_l / N_d #holdup correlation group

    ε_l_by_ψ = sqrt((0.0047 + 1123.32 * H + 729489.64 * H^2)/(1 + 1097.1566 * H + 722153.97 * H^2))

    B = N_gv * N_l^0.38 / N_d^2.14
    ψ = (1.0886 - 69.9473*B + 2334.3497*B^2 - 12896.683*B^3)/(1 - 53.4401*B + 1517.9369*B^2 - 8419.8115*B^3) #economides et al 235

    return ψ * ε_l_by_ψ
end


"""
HagedornAndBrownPressureDrop(pressure_est, id, v_sl, v_sg, ρ_l, ρ_g, μ_l, μ_g, σ_l, id_ft, λ_l, md, tvd, frictionfactor::Function, uphill_flow, roughness)

Helper function for H&B -- compute H&B pressure drop when bubble flow criteria are not met.
"""
function HagedornAndBrownPressureDrop(pressure_est, id, v_sl, v_sg, ρ_l, ρ_g, μ_l, μ_g, σ_l, id_ft, λ_l, md, tvd, frictionfactor::Function, uphill_flow,  roughness)

    ε_l = HagedornAndBrownLiquidHoldup(pressure_est, id, v_sl, v_sg, ρ_l, μ_l, σ_l)

    if uphill_flow
        ε_l = max(ε_l, λ_l) #correction to original: for uphill flow, true holdup must by definition be >= no-slip holdup
    end

    ρ_m = ρ_l * ε_l + ρ_g * (1 - ε_l) #mixture density in lb/ft³
    massflow = π*(id_ft/2)^2 * (v_sl * ρ_l + v_sg * ρ_g) * 86400 #86400 s/day

    #%% friction factor:
    μ_m = μ_l^ε_l * μ_g^(1-ε_l)
    N_Re = 2.2e-2 * massflow / (id_ft * μ_m)
    fric = frictionfactor(N_Re, id, roughness)/4 #corrected friction factor

    #%% core calculation:
    dpdl_el = 1/144.0 * ρ_m #elevation component
    friction_effect = uphill_flow ? 1 : -1 #note that friction MUST act against the direction of flow
    dpdl_f = friction_effect * 1/144.0 * fric * massflow^2 / (7.413e10 * id_ft^5 *ρ_m) #frictional component: see Economides et al
    #ρ_ns = λ_l * ρ_l + λ_g * ρ_g
    #dpdl_f = 1.294e-3 * fric * ρ_ns^2 * v_m^2 / (ρ_m * id) #Takacs -- takes normal friction factor
    #dpdl_kinetic = 2.16e-4 * ρ_m * v_m * (dvm_dh) #neglected except with high mass flow rates

    dp_dl = dpdl_el * tvd + dpdl_f * md #+ dpdl_kinetic * md

    return dp_dl
end


"""
GriffithWallisPressureDrop(v_sl, v_sg, v_m, ρ_l, ρ_g, μ_l, id_ft, md, tvd, frictionfactor::Function, uphill_flow, roughness)

Helper function for H&B correlation -- compute Griffith pressure drop for bubble flow regime.
"""
function GriffithWallisPressureDrop(v_sl, v_sg, v_m, ρ_l, ρ_g, μ_l, id_ft, md, tvd, frictionfactor::Function, uphill_flow,  roughness)
    v_s = 0.8 #assumed slip velocity of 0.8 ft/s -- probably assumes gas in oil bubbles with no water cut or vice versa?
    ε_l = 1 - 0.5 * (1 + v_m / v_s - sqrt((1 + v_m/v_s)^2 - 4*v_sg/v_s))

    if uphill_flow
        ε_l = max(ε_l, λ_l) #correction to original: for uphill flow, true holdup must by definition be >= no-slip holdup
    end

    ρ_m = ρ_l * ε_l + ρ_g * (1 - ε_l) #mixture density in lb/ft³
    massflow = π*(id_ft/2)^2 * (v_sl * ρ_l + v_sg * ρ_g) * 86400 #86400 s/day

    N_Re = 2.2e-2 * massflow / (id_ft * μ_l) #Reynolds number
    fric = frictionfactor(N_Re, id, roughness)/4 #corrected friction factor

    dpdl_el = 1/144.0 * ρ_m #elevation component
    friction_effect = uphill_flow ? 1 : -1 #note that friction MUST act against the direction of flow
    dpdl_f = friction_effect* 1/144.0 * fric * massflow^2 / (7.413e10 * id_ft^5 * ρ_l * ε_l^2) #frictional component
    dpdl = dpdl_el * tvd + dpdl_f * md

    return dpdl
end




"""
HagedornAndBrown(<arguments>)

Calculates pressure drop for a single pipe segment using the Hagedorn & Brown 1965 method (with recent modifications), with optional Griffith and Wallis bubble flow corrections.

Returns a ΔP in psi.

Doesn't account for oil/water phase slip, and does not incorporate flow regimes distinctions outside of in/out of bubble flow. Originally developed for vertical wells.

As of release v0.9, assumes **outlet-defined** models only, i.e. top-down from wellhead; thus, uphill flow corresponds to producers and downhill flow to injectors.

For more information, see *Petroleum Production Systems* by Economides et al., or the the Fekete [reference on pressure drops](http://www.fekete.com/san/webhelp/feketeharmony/harmony_webhelp/content/html_files/reference_material/Calculations_and_Correlations/Pressure_Loss_Calculations.htm).


# Arguments
All arguments take U.S. field units.

- `md`: measured depth of the pipe segment, feet
- `tvd`: true vertical depth, feet
- `inclination`: inclination from vertical, degrees (e.g. vertical => 0)
- `id`: inner diameter of the pipe segment, inches
- `v_sl`: superficial liquid mixture velocity, ft/s
- `v_sg`: superficial gas velocity, ft/s
- `ρ_l`: liquid mixture density, lb/ft³
- `ρ_g`: gas density, lb/ft³
- `σ_l`: liquid/gas interfacial tension, centipoise
- `μ_l`: liquid mixture dynamic viscosity
- `μ_g`: gas dynamic viscosity
- `roughness`: pipe roughness, inches
- `pressure_est`: estimated average pressure of the pipe segment (needed to determine the kinetic effects component of the pressure drop)
- `frictionfactor::Function = SerghideFrictionFactor`: function used to determine the Moody friction factor
- `uphill_flow = true`: indicates uphill or downhill flow. It is assumed that the start of the 1D segment is an outlet and not an inlet
- `GriffithWallisCorrection = true`: indicates whether the Griffith and Wallis 1961 corrections should be applied to prevent overprediction of liquid holdup.
"""
function HagedornAndBrown(md, tvd, inclination, id,
                            v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g,
                            roughness, pressure_est, frictionfactor::Function = SerghideFrictionFactor,
                            uphill_flow = true, GriffithWallisCorrection = true)

    id_ft = id/12

    v_m = v_sl + v_sg

    #%% holdup:
    λ_l = v_sl / v_m
    λ_g = 1 - λ_l

    if GriffithWallisCorrection
        L_B = max(1.071 - 0.2218 * v_m^2 / id, 0.13) #Griffith bubble flow boundary
        if λ_g <  L_B
            dpdl = GriffithWallisPressureDrop(v_sl, v_sg, v_m, ρ_l, ρ_g, μ_l, id_ft, md, tvd, frictionfactor, uphill_flow, roughness)
        else
            dpdl = HagedornAndBrownPressureDrop(pressure_est, id, v_sl, v_sg, ρ_l, ρ_g, μ_l, μ_g, σ_l, id_ft, λ_l, md, tvd, frictionfactor, uphill_flow, roughness)
        end
    else #no correction
        dpdl = HagedornAndBrownPressureDrop(pressure_est, id, v_sl, v_sg, ρ_l, ρ_g, μ_l, μ_g, σ_l, id_ft, λ_l, md, tvd, frictionfactor, uphill_flow, roughness)
    end

    return dpdl
end #Hagedorn & Brown
