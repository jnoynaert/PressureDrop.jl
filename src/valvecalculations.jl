using PrettyTables

#TODO: docs
"""
Indicate orifice valves with an R-value and PTRO of 0.
"""
struct GasliftValves

    md::Array{Float64,1}
    PTRO::Array{Float64,1}
    R::Array{Float64,1}
    port::Array{Int8,1}

    function GasliftValves(md::Array, PTRO::Array, R::Array, port::Array)

        ports = try
            convert(Array{Int8,1}, port)
        catch
            throw(ArgumentError("Specify port sizes in integer 64ths inches, e.g. 16 for a quarter-inch port."))
        end

        if any(R .> 1) || any(R .< 0)
            throw(ArgumentError("R-values are the area ratio of the port to the bellows and must be in [0, 1]."))
        elseif any(R .> 0.2)
            @info "Large R-value(s) entered--validate valve entry data."
        end

        new(convert(Array{Float64,1}, md), convert(Array{Float64,1}, PTRO), convert(Array{Float64,1}, R), ports)
    end
end


"""
ThornhillCraver_gaspassage_simplified(P_td, P_cd, T_cd, portsize_64ths)

Thornhill-Craver gas throughput for square-edged orifice (optimistic since it assumes a fully open valve where the stem does not interfere with flow).

This simplified version assumes gas specific gravity at 0.65.

# Arguments
- `P_td`: tubing pressure at depth, **psia**
- `P_cd`: casing pressure at depth, psia
- `T_cd`: gas (casing fluid) temperature at depth, °F
- `portsize_64ths`: valve port size in 64ths inch
"""
function ThornhillCraver_gaspassage_simplified(P_td, P_cd, T_cd, portsize_64ths)

    if P_td > P_cd
        return 0
    end

    A_p = π * (portsize_64ths/128)^2 #port area, in²

    R_du = max(P_td / P_cd, 0.553) #critical flow condition check

    return 2946 * A_p * P_cd * sqrt(R_du^1.587 - R_du^1.794) / sqrt(T_cd + 459.67)
end


"""
ThornhillCraver_gaspassage(<args>)

Thornhill-Craver gas throughput for square-edged orifice.

See section 8.1 of *Fundamentals of Gas Lift Engineering* by Ali Hernandez, as well as published work by Ken Decker, for an in-depth discussion of the application of Thornhill-Craver to gas lift valve performance.

In general, T-C will be optimistic, but should be expected to have an effective error of up to +/- 30%.

# Arguments
- `P_td`: tubing pressure at depth, **psia**
- `P_cd`: casing pressure at depth, psia
- `T_cd`: gas (casing fluid) temperature at depth, °F
- `portsize_64ths`: valve port size in 64ths inch
- `sg_gas`: gas specific gravity relative to air
- `molFracCO2 = 0.0`, `molFracH2S = 0.0`: produced gas fractions of hydrogen sulfide and CO2, [0,1]
- `C_d = 0.827`: discharge coefficient--uses 0.827 by defaul to match original T-C
- `Z_correlation::Function = KareemEtAlZFactor`: natural gas compressibility/Z-factor correlation to use
- `pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure`: psuedocritical pressure function to use
- `pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp`: pseudocritical temperature function to use
"""
function ThornhillCraver_gaspassage(P_td, P_cd, T_cd, portsize_64ths,
                                    sg_gas, molFracCO2 = 0, molFracH2S = 0, C_d = 0.827, Z_correlation::Function = KareemEtAlZFactor, pseudocrit_pressure_correlation::Function = HankinsonWithWichertPseudoCriticalPressure, pseudocrit_temp_correlation::Function = HankinsonWithWichertPseudoCriticalTemp)

    if P_td > P_cd
        return 0
    end

    k = 1.27 #assumed value for Cp_Cv; if a more precise solution is needed later, see SPE 150808 "Specific Heat Capacity of Natural Gas; Expressed as a Function of ItsSpecific gravity and Temperature" by Kareem Lateef et al

    P_pc = pseudocrit_pressure_correlation(sg_gas, molFracCO2, molFracH2S)
    _, T_pc, _ = pseudocrit_temp_correlation(sg_gas, molFracCO2, molFracH2S)
    Z = Z_correlation(P_pc, T_pc, P_cd, T_cd) #compressibility factor at upstream pressure

    T_gd = T_cd + 459.67
    A_p = π * (portsize_64ths/128)^2 #port area, in²
    F_cf = (2/(k+1))^(k/(k-1)) #critical flow ratio
    F_du = max(P_td/P_cd, F_cf) #critical flow condition check

    return C_d * A_p * 155.5 *  P_cd * sqrt(2 * 32.16 * k / (k-1) * (F_du^(2/k) - F_du^((k+1)/k))) / sqrt(Z * sg_gas * T_gd)
end


"""
Sage and Lacy experimental Z-factor correlation.

Takes pressure in psia and temperature in °F.
"""
function SageAndLacy_nitrogen_Zfactor(p, T)
    b = (1.207e-7 * T^3 - 1.302e-4 * T^2 + 5.122e-2 * T - 4.781) * 1e-5
    c = (-2.461e-8 * T^3 + 2.640e-5 * T^2 - 1.058e-2 * T + 1.880) * 1e-8

    return 1 + b * p + c * p^2
end


"""
domepressure_downhole(p_d_set, T_v, error_tolerance = 0.1, p_d_est = p_d_set/0.9, T_set = 60, Zfactor::Function = SageAndLacy_nitrogen_Zfactor)

Iteratively calculates the dome pressure of the valve downhole using gas equation of state, assuming that the change in dome volume is negligible.

# Arguments

- `PTRO`: test rack opening pressure of valve **in psig**
- `R`: R-value of the valve (nominally, area of port divided by area of bellows/dome, but adjusted for lapped seats, etc)
- `T_v`: valve temperature at depth
- `error_tolerance = 0.1`: error tolerance in psi
- `delta_est`: initial estimate for downhole dome pressure as a percentage of surface set pressure
- `T_set = 60`: setting temperature in °F
- `Zfactor::Function = SageAndLacy_nitrogen_Zfactor`: Z-factor function parameterized by target conditions as pressure in psia and temperature in °F
"""
function domepressure_downhole(PTRO, R, T_v, error_tolerance = 0.1, delta_est = 1.1, T_set = 60, Zfactor::Function = SageAndLacy_nitrogen_Zfactor)

    p_d_set = (PTRO + 14.65) * (1 - R)
    p_d_est = p_d_set * delta_est
    p_d = Zfactor(p_d_est, T_v) / Zfactor(p_d_set, T_set) * (T_v + 459.67) * p_d_set / 520

    while abs(p_d - p_d_est) > error_tolerance
        p_d_est = p_d
        p_d = Zfactor(p_d_est, T_v) / Zfactor(p_d_set, T_set) * (T_v + 459.67) * p_d_set / 520
    end

    return p_d
end


#TODO: clean doc
"""

**CAUTION**: All inputs are in psia, but PSO/PSCs are returned in psig for ease of interpretation.

All forms are derived from the force balance for opening, P_t * A_p + P_c * (A_d - A_p) ≥ P_d * A_d
    and the force balance for closing, P_c ≤ P_d.

Note that:
- the valve closing pressures given are a theoretical minimum (casing pressure is assumed to act on the entire area of the bellows/dome until closing).
- valve opening and closing pressures are recalculated from PTRO and current conditions, rather than the other way around common during design.
- casing temperature is adjusted from tubing temperature if only a tubing temperature profile is provided. **To force the use of identical temperature profiles, pass the tubing temperature twice.**
"""
function valve_calcs(valves::GasliftValves, well::Wellbore, sg_gas, tubing_pressures::Array{T, 1} where T <: Real, casing_pressures::Array{T, 1} where T <: Real,
                    tubing_temps::Array{T, 1} where T <: Real, casing_temps::Array{T, 1} where T <: Real = tubing_temps .* 0.85,
                    dp_min = 100, one_inch_coefficient = 0.76, one_pt_five_inch_coefficient = 0.6)

    interp_values = interpolate_all(well, [tubing_pressures, casing_pressures, tubing_temps, casing_temps, well.tvd], valves.md)

    P_td = interp_values[:,1] #tubing pressure at depth
    P_cd = interp_values[:,2] #casing pressure at depth

    T_td = interp_values[:,3] #temp of tubing fluid at depth
    T_cd = interp_values[:,4] #temp of casing fluid at depth

    valve_temps = 0.7 .* T_cd .+ 0.3 .* T_td # Faustinelli, J.G. 1997

    PVC = P_bd = domepressure_downhole.(valves.PTRO, valves.R, valve_temps) #dome/bellows pressure at depth; equals valve closing pressure at depth

    PVO = (P_bd .- (P_td .* valves.R)) ./ (1 .- valves.R) #valve opening pressure at depth

    csg_diffs = P_cd .- casing_pressures[1] #casing differential to depth
    PSC = PVC .- csg_diffs
    PSO = PVO .- csg_diffs

    T_C = ThornhillCraver_gaspassage.(P_td, P_cd, T_cd, valves.port, sg_gas)

    GLV_numbers = length(valves.md):-1:1

    PPEF = valves.R ./ (1 .- valves.R) .* 100 #production pressure effect factor

    # returns PSOs in psig to avoid confusion
    valvetable = hcat(GLV_numbers, valves.md, interp_values[:,5], PSO .- 14.65, PSC .- 14.65, valves.port, valves.R, PPEF, valves.PTRO,
                P_td, P_cd, PVO, PVC, T_td, T_cd, T_C, T_C * one_inch_coefficient, T_C * one_pt_five_inch_coefficient)

    #lowest valve where PVO <= CP && PVC > CP && R-value != 0
    valve_ids = collect(1:length(valves.md))
    active_valve_row = findlast(i -> P_cd[i] - P_td[i] >= dp_min && valves.R[i] > 0 && PVO[i] <= P_cd[i] && PVC[i] > P_cd[i], valve_ids) #only check non-orifice valves
    if active_valve_row === nothing
        if P_cd[end] - P_td[end] >= dp_min #operating on orifice valve
            active_valve_row = length(valves.md)
        else
            active_valve_row = 1 #nominal lockout, but assume you are injecting on top valve
        end
    end
    injection_depth = valves.md[active_valve_row]

    return valvetable, injection_depth
end


#TODO: add doc
"""
"""
function valve_table(valvedata, injection_depth)
    header = ["GLV" "MD" "TVD" "PSO" "PSC" "Port" "R" "PPEF" "PTRO" "TP" "CP" "PVO" "PVC" "T_td" "T_cd" "Q_o" "Q_1.5" "Q_1";
                "" "ft" "ft" "psig" "psig" "64ths" "" "%" "psig" "psia" "psia" "psia" "psia" "°F" "°F" "mcf/d" "mcf/d" "mcf/d"]

    operating_valve = Highlighter((data, i, j) -> (data[i,2] == injection_depth), crayon"bg:dark_gray white bold")

    pretty_table(valvedata, header, unicode_rounded; header_crayon = crayon"yellow bold",
                 formatter = ft_printf("%.0f", [1:6;8:18]), highlighters = (operating_valve,))
end


"""
estimate_valve_Rvalue(port_size, valve_size, lapped_seat = true)

Estimates an R-value for your valve (not recommended) using sensible defaults, if you do not have a manufacturer-provided value specific to the valve (recommended).

Takes port size as a diameter in 64ths inches, a valve size in inches {1.0, 1.5}, and an indication of whether the seat is lapped as {true, false}.
"""
function estimate_valve_Rvalue(port_size, valve_size, lapped_seat = true)

    if !(valve_size ∈ (1.0, 1.5))
        throw(ArgumentError("Must specify a 1.0\" or 1.5\" valve."))
    end

    port_diameter = lapped_seat ? (port_size/64 + 0.006) : port_size
    port_area = π * (port_diameter/2)^2
    bellows_area = valve_size == 1.0 ? 0.31 : 0.77

    return port_area/bellows_area
end
