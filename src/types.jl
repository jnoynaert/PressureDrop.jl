"""
struct GasliftValves: a type to define a string of gas lift valves for valve & pressure calculations.

Constructor: `GasliftValves(md::Array, PTRO::Array, R::Array, port::Array)`

Port sizes must be in integer increments of 64ths inches.

Indicate orifice valves with an R-value and PTRO of 0.
"""
struct GasliftValves

    md::Array{Float64,1}
    PTRO::Array{Float64,1}
    R::Array{Float64,1}
    port::Array{Int8,1}

    function GasliftValves(md::Array{T} where T <: Real, PTRO::Array{T} where T <: Real, R::Array{T} where T <: AbstractFloat, port::Array{T} where T <: Union{Real, Int})

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


#printing for gas lift valves
Base.show(io::IO, valves::GasliftValves) = print(io, "Valve design with $(length(valves.md)) valves and bottom valve at $(valves.md[end])' MD.")


"""
struct Wellbore: type to define a flow path as an input for pressure drop calculations

See `read_survey` for helper method to create a Wellbore object from deviation survey files.

# Fields
- `md::Array{Float64, 1}`: measured depth for each segment in feet
- `inc::Array{Float64, 1}`: inclination from vertical for each segment in degrees, e.g. true vertical = 0°
- `tvd::Array{Float64, 1}`: true vertical depth for each segment in feet
- `id::Array{Float64, 1}`: inner diameter for each pip segment in inches

# Constructors
By default, negative depths are disallowed, and a 0 MD / 0 TVD point is added if not present, to allow graceful handling of outlet pressure definitions.
To bypass both the error checking and convenience feature, pass `true` as the final argument to the constructor.

`Wellbore(md, inc, tvd, id::Array{Float64, 1}, allow_negatives = false)`: defines a new Wellbore object from a survey with inner diameter defined for each segment. Lengths of each input array must be equal.

`Wellbore(md, inc, tvd, id::Float64, allow_negatives = false)`: defines a new Wellbore object with a uniform ID along the entire flow path.

`Wellbore(md, inc, tvd, id, valves::GasliftValves, allow_negatives = false)`: defines a new Wellbore object and adds interpolated survey points for each gas lift valve.
"""
struct Wellbore

    md::Array{Float64, 1}
    inc::Array{Float64, 1}
    tvd::Array{Float64, 1}
    id::Array{Float64, 1}

    function Wellbore(md, inc, tvd, id::Array{Float64, 1}, allow_negatives::Bool = false)

        lens = length.([md, inc, tvd, id])

        if !( count(x -> x == lens[1], lens) == length(lens) )
            throw(DimensionMismatch("Mismatched number of wellbore elements used in wellbore constructor."))
        end

        if !allow_negatives
            if minimum(md) < 0 || minimum(tvd) < 0
                throw(ArgumentError("Survey contains negative measured or true vertical depths. Pass the `allow_negatives` constructor flag if this is intentional."))
            end

            #add the origin/outlet reference point if missing
            if !(md[1] == tvd[1] == 0)
                md = vcat(0, md)
                inc = vcat(0, inc)
                tvd = vcat(0, tvd)
                id = vcat(id[1], id)
            end
        end

        new(md, inc, tvd, id)
    end
end #struct Wellbore


#convenience constructor for uniform tubulars
Wellbore(md, inc, tvd, id::Float64, allow_negatives::Bool = false) = Wellbore(md, inc, tvd, repeat([id], inner = length(md)), allow_negatives)


#convenience constructors to add reference depths for valves so that they can be used as injection points
#TODO: add an error message to the injection point logic that makes it clear how to fix the problem if they don't use this constructor and you get an injection point that doesn't match the wellbore
function Wellbore(md, inc, tvd, id, valves::GasliftValves, allow_negatives::Bool = false)

    well = Wellbore(md, inc, tvd, id, allow_negatives)

    for v in 1:length(valves.md)

        upper_index = searchsortedlast(well.md, valves.md[v])

        if well.md[upper_index] != valves.md[v]
            lower_index = upper_index + 1 #also the target insertion position

            x1, x2 = well.md[upper_index], well.md[lower_index]
            for property in [well.inc, well.tvd, well.id]
                y1, y2 = property[upper_index], property[lower_index]
                interpolated_value = y1 + (y2 - y1)/(x2 - x1) * (valves.md[v] - x1)
                insert!(property, lower_index, interpolated_value)
            end

            insert!(well.md, lower_index, valves.md[v])
        end
    end

    return well
end


#Printing for Wellbore structs
Base.show(io::IO, well::Wellbore) = print(io,
    "Wellbore with $(length(well.md)) points.\n",
    "Ends at $(well.md[end])' MD / $(well.tvd[end])' TVD.\n",
    "Max inclination $(maximum(well.inc))°. Average ID $(round(sum(well.id)/length(well.id), digits = 3)) in.")


#model struct. ONLY applies to wrapper functions.
#TODO: docs
"""
Makes it easier to iterate well models
"""
mutable struct WellModel

    wellbore::Wellbore; roughness
    valves::Union{GasliftValves, Missing}
    temperatureprofile::Union{Array{T,1}, Missing} where T <: Real
    temperature_method; WHT; geothermal_gradient; BHT; casing_temp_factor
    pressurecorrelation::Function; outlet_referenced::Bool
    WHP; CHP; dp_est; dp_est_inj; error_tolerance; error_tolerance_inj
    q_o; q_w; GLR
    injection_point; naturalGLR
    APIoil; sg_water; sg_gas; sg_gas_inj
    molFracCO2; molFracH2S; molFracCO2_inj; molFracH2S_inj
    pseudocrit_pressure_correlation::Function
    pseudocrit_temp_correlation::Function
    Z_correlation::Function
    gas_viscosity_correlation::Function
    solutionGORcorrelation::Function
    oilVolumeFactor_correlation::Function
    waterVolumeFactor_correlation::Function
    dead_oil_viscosity_correlation::Function
    live_oil_viscosity_correlation::Function
    frictionfactor::Function

    function WellModel(;wellbore, roughness, valves = missing, temperatureprofile = missing,
                        temperature_method = "linear", WHT = missing, geothermal_gradient = missing, BHT = missing, casing_temp_factor = 0.85,
                        pressurecorrelation = BeggsAndBrill, outlet_referenced = true,
                        WHP, CHP = missing, dp_est, dp_est_inj = 0.1 * dp_est, error_tolerance = 0.1, error_tolerance_inj = 0.05,
                        q_o, q_w, GLR, injection_point = missing, naturalGLR = missing,
                        APIoil, sg_water, sg_gas, sg_gas_inj = sg_gas,
                        molFracCO2 = 0.0, molFracH2S = 0.0, molFracCO2_inj = molFracCO2, molFracH2S_inj = molFracH2S,
                        pseudocrit_pressure_correlation = HankinsonWithWichertPseudoCriticalPressure,
                        pseudocrit_temp_correlation = HankinsonWithWichertPseudoCriticalTemp,
                        Z_correlation = KareemEtAlZFactor, gas_viscosity_correlation = LeeGasViscosity,
                        solutionGORcorrelation = StandingSolutionGOR, oilVolumeFactor_correlation = StandingOilVolumeFactor,
                        waterVolumeFactor_correlation = GouldWaterVolumeFactor,
                        dead_oil_viscosity_correlation = GlasoDeadOilViscosity, live_oil_viscosity_correlation = ChewAndConnallySaturatedOilViscosity,
                        frictionfactor = SerghideFrictionFactor)

        new(wellbore, roughness, valves, temperatureprofile, temperature_method, WHT, geothermal_gradient, BHT, casing_temp_factor,
            pressurecorrelation, outlet_referenced, WHP, CHP, dp_est, dp_est_inj, error_tolerance, error_tolerance_inj,
            q_o, q_w, GLR, injection_point, naturalGLR, APIoil, sg_water, sg_gas, sg_gas_inj, molFracCO2, molFracH2S, molFracCO2_inj, molFracH2S_inj,
            pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation, gas_viscosity_correlation, solutionGORcorrelation,
            oilVolumeFactor_correlation, waterVolumeFactor_correlation, dead_oil_viscosity_correlation, live_oil_viscosity_correlation, frictionfactor)
    end

end

#Printing for model structs
function Base.show(io::IO, model::WellModel)

    fields = fieldnames(WellModel)
    values = map(f -> getfield(model, f), fields) |> list -> map(x -> !(x isa Array) ? string(x) : "$(length(x)) points from $(maximum(x)) to $(minimum(x)).", list)
    msg = string.(fields) .* " : " .* values

    println(io, "Well model: ")
    for item in msg
        println(io, item)
    end
end
