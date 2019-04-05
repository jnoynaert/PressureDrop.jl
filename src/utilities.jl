"""
Assumes MD, INC, TVD, ID order for columns.

Skip header of n lines with skiplines = n.
"""
function read_survey(;path::String, delim::Char = ',', skiplines::Int64 = 1, maxdepth::Union{Bool, Real} = false, id_included = false, id = 2.441)

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

    return id_included ? Wellbore(output[:,1], output[:,2], output[:,3], output[:,4]) :
        Wellbore(output[:,1], output[:,2], output[:,3], id)
end
