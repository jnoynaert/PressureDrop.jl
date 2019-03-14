# Pressure correlations for PressureDrop package.

#TODO: ensure appropriate kwarg handling so that all correlations utilize a common interface.


#%% Helper functions

"""
takacs 52

Note that this does not account for slip between liquid phases.
"""
function liquidvelocity_superficial(q_o, q_w, id, B_o, B_w)
    A = π * (id/2.0)^2

    if q_o > 0
        WOR = q_w / q_o
        return 6e-5 * (q_o + q_w) / A * (B_o/(1 + WOR) + B_w * WOR / (1 + WOR))
    else #100% WC
        return q_w * B_w / A
    end
end


"""
takacs 52

"""
function gasvelocity_superficial(q_o, q_w, GLR, R_s, id, B_g)
    A = π * (id/2.0)^2

    if q_o > 0
        WOR = q_w / q_o
        return 1.16e-5 * (q_o + q_w) / A * (GLR - R_s /(1 + WOR)) * B_g
    else #100% WC
        return 1.16e-5 * q_w * (GLR - R_s) * B_g / A
    end
end

# mixture velocity: just v_sg + v_sl


"""
Weighted average for mixture properties.

Does not account for oil slip, mixing effects, fluid expansion, Non-Newtonian behavior of emulsions, etc.
"""
function mixture_properties_simple(q_o, q_w, p_o, p_w)

    return (q_o * p_o + q_w * p_w) / (q_o + q_w)
end


"""
k is epsilon/d, where epsilon = 0.0006ish for new and 0.009ish for used
and epsilon is the absolute roughness in inches
and k is the pipe ID in inches.

may only apply for turbulent flow? otherwise will overpredict the loss

Takacs p30
"""
function ChenFrictionFactor(N_Re, id, roughness = 0.01) #TODO: verify this and set up a test

    k = roughness/id

    A = k^1.1098 / 2.8257 + (7.149 / N_Re)^0.8981
    x = -2 * log10(k / 3.7065 - 5.0452 / N_Re * log10(A))

    return 1/x^2
end
#TODO: handle laminar flow


#%% Beggs and Brill

"""
Beggs and Brill flow pattern as a string ∈ {"segregated", "transition", "distributed", "intermittent"}.

Takes no-slip holdup (λ_l) and mixture Froude number (N_Fr).

Beggs and Brill. Takacs p87.
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
#TODO: see flow pattern definitions at https://wiki.pengtools.com/index.php?title=Beggs_and_Brill_correlation#cite_note-BB1991-2, which
# appear to be more robust.



"""
Helper function for Beggs and Brill. Returns adjusted liquid holdup, ε_l(α).

Optional Payne et al correction applied to output.

Takes flow pattern (string ∈ {"segregated", "intermittent", "distributed"}), no-slip holdup (λ_l), Froude number (N_Fr),
liquid velocity number (N_lv), angle from horizontal (α, radians), uphill flow (boolean).

Ref Takacs 88.
"""
function BeggsAndBrillAdjustedLiquidHoldup(flowpattern, λ_l, N_Fr, N_lv, α, inclination, uphill_flow, PayneCorrection = true) #TODO: add a test

    if PayneCorrection && uphill_flow
        correctionfactor = 0.924
    elseif PayneCorrection
        correctionfactor = 0.685
    else
        correctionfactor = 1.0
    end

    BB_coefficients = Dict("segregated" => (a = 0.980, b= 0.4846, c = 0.0868, e = 0.011, f = -3.7680, g = 3.5390, h = -1.6140),
                            "intermittent" => (a = 0.845, b = 0.5351, c = 0.0173, e = 2.960, f = 0.3050, g = -0.4473, h = 0.0978),
                            "distributed" => (a = 1.065, b = 0.5824, c = 0.0609),
                            "downhill" => (e = 4.700, f = -0.3692, g = 0.1244, h = -0.5056) )

    a = BB_coefficients[flowpattern][:a]
    b = BB_coefficients[flowpattern][:b]
    c = BB_coefficients[flowpattern][:c]

    ε_l_horizontal = a * λ_l^b / N_Fr^c #liquid holdup assuming horizontal (α = 0 rad)

    #TODO: add verification (compare horizontal to noslip, see steps on Takacs 90)

    if α ≈ 0 #horizontal flow
        return ε_l_horizontal
    else #inclined or vertical flow
        if uphill_flow
            if flowpattern == "distributed"
                ψ = 1.0
            else
                e = BB_coefficients[flowpattern][:e]
                f = BB_coefficients[flowpattern][:f]
                g = BB_coefficients[flowpattern][:g]
                h = BB_coefficients[flowpattern][:h]

                C = max( (1 - λ_l) * log(e * λ_l^f * N_lv^g * N_Fr^h), 0)

                if inclination ≈ 0 #vertical flow
                    ψ = 1 + 0.3 * C
                else
                    ψ = 1 + C * (sin(1.8*α) - (1/3) * sin(1.8*α)^3)
                end
            end
        else #downhill flow
            e = BB_coefficients["downhill"][:e]
            f = BB_coefficients["downhill"][:f]
            g = BB_coefficients["downhill"][:g]
            h = BB_coefficients["downhill"][:h]

            C = max( (1 - λ_l) * log(e * λ_l^f * N_vl^g * N_Fr^h), 0)

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
Documentation here.

Doesn't account for oil/water phase slip.

Currently assumes *outlet-defined* models only, i.e. top-down from wellhead (easy because they always converge to an inlet/BHP); thus, uphill flow corresponds to producers and downhill flow to injectors.

Returns a ΔP in psi.

http://www.fekete.com/san/webhelp/feketeharmony/harmony_webhelp/content/html_files/reference_material/Calculations_and_Correlations/Pressure_Loss_Calculations.htm
for additional ref
"""
function BeggsAndBrill( md, tvd, inclination, id,
                        v_sl, v_sg, ρ_l, ρ_g, σ_l, μ_l, μ_g,
                        roughness, pressure_est,
                        uphill_flow = true, PayneCorrection = true )

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

    #%% friction factor:
    y = λ_l / ε_l_adj^2
    if 1.0 < y < 1.2
        s = ln(2.2y - 1.2) #handle the discontinuity
    else
        ln_y = log(y)
        s = ln_y / (-0.0523 + 3.182 * ln_y - 0.872 * ln_y^2 + 0.01853 * ln_y^4)
    end

    fbyfn = ℯ^s #f/fₙ

    ρ_ns = ρ_l * λ_l + ρ_g * (1-λ_l) #no-slip density
    μ_ns = μ_l * λ_l + μ_g * (1-λ_l) #no-slip friction in centipoise
    N_Re = 124 * ρ_ns * v_m * id / μ_ns #Reynolds number
    f_n = ChenFrictionFactor(N_Re, id, roughness)
    fric = f_n * fbyfn #friction factor


    #%% core calculation:
    ρ_m = ρ_l * ε_l_adj + ρ_g * (1 - ε_l_adj) #mixture density in lb/ft³

    dpdl_el = (1/144.0) * ρ_m * sin(α) #elevation component
    friction_effect = uphill_flow ? 1 : -1 #note that friction MUST act against the direction of flow
    dpdl_f = friction_effect * 1.294e-3 * fric * (ρ_ns * v_m^2) / id #frictional component
    E_k = 2.16e-4 * fric * (v_m * v_sg * ρ_ns) / pressure_est #kinetic effects; typically negligible

    dp_dl = (dpdl_el * tvd + dpdl_f * md) / (1 - friction_effect*E_k) #assumes friction and kinetic effects both increase pressure in the same 1D direction

    return dp_dl
end #TODO: add tests
