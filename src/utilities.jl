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

`path::String`: absolute or relative path to survey file
`delim::Char = ','`: file delimiter
`skiplines::Int64 = 1`: number of lines to skip before survey data starts; assumes a 1-line header by default
`maxdepth::Union{Bool, Real} = false`: If set to a real number, drop survey data past a certain measured depth. If false, keep the entire survey.
`id_included::Bool = false`: whether the survey segment ID is stored as a fourth column. This is the easiest option to include tapered strings.
`id::Real = 2.441`: the diameter to assume for the entire wellbore length, if the ID is not included in the survey file.
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
