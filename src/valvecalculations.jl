using PrettyTables


# maintain strict typing to protect type inference when reading configurations from files
struct GasliftValves

    md::Array{Float64,1}
    PTRO::Array{Float64,1}
    R::Array{Float64,1}
    port::Array{Int8,1}

    #TODO: constructor that casts inputs correctly
end


"""
"""
function ThornhillCraver_Q()
    #TODO: placeholder
    return 1
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

- `PTRO`: test rack opening pressure of valve
- `R`: R-value of the valve (nominally, area of port divided by area of bellows/dome, but adjusted for lapped seats, etc)
- `T_v`: valve temperature at depth
- `error_tolerance = 0.1`: error tolerance in psi
- `p_d_est = p_d_set/0.9`: initial estimate for downhole dome pressure
- `T_set = 60`: setting temperature in °F
- `Zfactor::Function = SageAndLacy_nitrogen_Zfactor`: Z-factor function parameterized by target conditions as pressure in psia and temperature in °F
"""
function domepressure_downhole(PTRO, R, T_v, error_tolerance = 0.1, p_d_est = p_d_set/0.9, T_set = 60, Zfactor::Function = SageAndLacy_nitrogen_Zfactor)

    p_d_set = PTRO * (1 - R)
    p_d = Zfactor(p_d_est, T_v) / Zfactor(p_d_set, T_set) * (T_v + 459.67) * p_d_set / 520

    while abs(p_d - p_d_est) > error_tolerance
        p_d_est = p_d
        p_d = Zfactor(p_d_est, T_v) / Zfactor(p_d_set, T_set) * (T_v + 459.67) * p_d_set / 520
    end

    return p_d
end


"""

**CAUTION**: All inputs are in psia, but PSO/PSCs are returned in psig for ease of interpretation.

All forms are derived from the force balance for opening, P_t * A_p + P_c * (A_d - A_p) ≥ P_d * A_d
    and the force balance for closing, P_c ≤ P_d.

Note that:
- the valve closing pressures given are a theoretical minimum (casing pressure is assumed to act on the entire area of the bellows/dome until closing).
- valve opening and closing pressures are recalculated from PTRO and current conditions, rather than the other way around common during design.
- casing temperature is adjusted from tubing temperature if only a tubing temperature profile is provided. **To force the use of identical temperature profiles, pass the tubing temperature twice.**
"""
function valve_calcs(valves::GasliftValves, well::Wellbore, tubing_pressures::Array{Real, 1}, casing_pressures::Array{Real, 1},
                    tubing_temps::Array{Real, 1}, casing_temps::Array{Real, 1} = tubing_temps .* 0.85,
                    one_inch_coefficient = 0.8, one_pt_five_inch_coefficient = 0.6)

    interp_values = interpolate_all(well, [tubing_pressures, casing_pressures, tubing_temps, casing_temps, well.tvd], valves.md)

    P_td = interp_values[:,1] #tubing pressure at depth
    P_cd = interp_values[:,2] #casing pressure at depth

    T_td = interp_values[:,3] #temp of tubing fluid at depth
    T_cd = interp_values[:,4] #temp of casing fluid at depth

    valve_temps = 0.7 .* T_cd .+ 0.3 .* T_td # Faustinelli, J.G. 1997

    PVC = P_bd = domepressure_downhole.(valves.PTRO, valves.R, valve_temps) #dome/bellows pressure at depth; equals valve closing pressure at depth

    PVO = (P_bd .- (P_td .* valves.R)) ./ (1 .- valves.R) #valve opening pressure at depth

    csg_diffs = P_cd .- casing_pressures[1] #casing differential to depth
    PSC = PVC .+ csg_diffs
    PSO = PVO .+ csg_diffs

    T_C = [ThornhillCraver_Q() for valve in valves.md] #TODO: placeholder

    GLV_numbers = length(valves):-1:1

    # returns PSOs in psig to avoid confusion
    return hcat(GLV_numbers, valves.md, interp_values[:,5], PSO .- 14.65, PSC .- 14.65, valves.port, valves.R, valves.PTRO, P_td, P_cd, PVO, PVC, T_td, T_cd, )
end


#TODO: switch this so that valve_table wraps valve_calcs
"""
"""
function valve_table(valvedata)
    header = ["GLV" "MD" "TVD" "PSO" "PSC" "Port" "R" "PTRO" "TP" "CP" "PVO" "PVC" "T_td" "T_cd" "Q_i_max" "Q_i_1\"" "Q_i_1.5\"";
                "" "ft" "ft" "psig" "psig" "64ths in", "Ap/Ab", "°F", "psia", "psia", "psia", "psia", "°F", "°F", "mcf/d", "mcf/d", "mcf/d"]


    pretty_table(valvedata, header, unicode_rounded; header_crayon = crayon"yellow bold", formatter = ft_printf("%5f"))
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
