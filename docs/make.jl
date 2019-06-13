#push!(LOAD_PATH,"../src/")

using Documenter, PressureDrop

makedocs(sitename="PressureDrop.jl", format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")

#TODO: see how to set up automatic doc deployment on Travis CI https://github.com/JuliaDocs/Documenter.jl/blob/master/docs/src/man/hosting.md

deploydocs(
    rep = "github.com/jnoynaert/PressureDrop.jl.git"
)
