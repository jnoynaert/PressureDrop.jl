# utility functions for PressureDrop package

"""
read_survey(<named arguments>)

Reads in a wellbore deviation survey from a delimited file and returns a Wellbore object for use with pressure drop calculations.

Assumes the column order for the file is (MD, INC, TVD, <optional ID>), in U.S. field units, where:
    MD = measured depth, ft
    Inc = inclination, degrees from vertical
    TVD = true vertical depth, ft
    ID = pipe hydraulic diameter, inches

# Arguments

- `path::String`: absolute or relative path to survey file
- `delim::Char = ','`: file delimiter
- `skiplines::Int64 = 1`: number of lines to skip before survey data starts; assumes a 1-line header by default
- `maxdepth::Union{Bool, Real} = false`: If set to a real number, drop survey data past a certain measured depth. If false, keep the entire survey.
- `id_included::Bool = false`: whether the survey segment ID is stored as a fourth column. This is the easiest option to include tapered strings.
- `id::Real = 2.441`: the diameter to assume for the entire wellbore length, if the ID is not included in the survey file.
"""
function read_survey(;path::String, delim::Char = ',', skiplines::Int64 = 1, maxdepth::Union{Bool, Real} = false, id_included::Bool = false, id::Real = 2.441, allow_negatives::Bool = false)

    nlines = countlines(path) - skiplines
    ncols = id_included ? 4 : 3
    output = Array{Float64, 2}(undef, nlines, ncols)
    filestream = open(path, "r")

    try
        for skip in 1:skiplines
            readline(filestream)
        end

        for (index, line) in enumerate(eachline(filestream))
            parsedline = parse.(Float64, split(line, delim, keepempty = false))
            output[index, :] = parsedline
        end
    finally
        close(filestream)
    end

    if maxdepth != false
        output = output[output[:,1] .<= maxdepth, :]
    end

    return id_included ? Wellbore(output[:,1], output[:,2], output[:,3], output[:,4], allow_negatives) :
        Wellbore(output[:,1], output[:,2], output[:,3], id, allow_negatives)
end


"""
read_valves(;path::String, delim::Char = ',', skiplines::Int64 = 1)

Expects a CSV with columns for [measured depth (ft)], [test rack opening pressure (psig)], [R-value (dimensionless)], [port size (diameter in 64ths inch)].

Indicate orifice valves with an R-value and PTRO of 0.

# Arguments
- `path::String`: absolute or relative path to survey file
- `delim::Char = ','`: file delimiter
- `skiplines::Int64 = 1`: number of lines to skip before survey data starts; assumes a 1-line header by default
"""
function read_valves(;path::String, delim::Char = ',', skiplines::Int64 = 1)

    nlines = countlines(path) - skiplines
    output = Array{Float64, 2}(undef, nlines, 4)
    filestream = open(path, "r")

    try
        for skip in 1:skiplines
            readline(filestream)
        end

        for (index, line) in enumerate(eachline(filestream))
            parsedline = parse.(Float64, split(line, delim, keepempty = false))
            output[index, :] = parsedline
        end
    finally
        close(filestream)
    end

    return GasliftValves(output[:,1], output[:,2], output[:,3], output[:,4]) #constructor will parse appropriately
end


"""
interpolate(ref_array, property::Array{Real,1}, point)

Interpolate between points without any bounds checking.
"""
function interpolate(ref_array::Array{T,1} where T <: Real, property::Array{T,1} where T <: Real, point)

    index_above = searchsortedlast(ref_array, point)
    index_below = index_above + 1
    x1 = ref_array[index_above]

    if x1 == point
        interpolated_value = property[index_above]
    else
        x2 = ref_array[index_below]
        y1, y2 = property[index_above], property[index_below]
        interpolated_value = y1 + (y2 - y1)/(x2 - x1) * (point - x1)
    end

    return interpolated_value
end


"""
interpolate_all(well::Wellbore, properties::Array{Array{Real,1},1}, points::Array{Real,1})

Interpolate multiple points for multiple properties with bounds checking.
"""
function interpolate_all(well::Wellbore, properties::Array{Array{T,1},1} where T <: Real, points::Array{T,1} where T <: Real)

    if !(all(length.(properties) .== length(well.md)))
        throw(DimensionMismatch("Property array lengths and number of wellbore survey points must match."))
    end

    if !(all(points .<= well.md[end]) && all(points .>= well.md[1]))
        throw(DimensionMismatch("Interpolation points cannot be outside wellbore measured depth."))
    end

    results = Array{Float64, 2}(undef, length(points), length(properties))
    for (index, property) in enumerate(properties)
        results[:, index] = map(pt -> interpolate(well.md, property, pt), points)
    end

    return results
end
