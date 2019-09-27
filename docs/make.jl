using Documenter, PressureDrop

isCI = get(ENV, "CI", nothing) == "true" #Travis populates this env variable by default

makedocs(
    clean = true,
    pages = [   "Overview" => "index.md",
                "Core functions" => "core.md", #includes types
                "Plotting" => "plotting.md",
                "Pressure & temperature correlations" => "correlations.md",
                "PVT properties" => "pvt.md",
                "Valve calculations" => "valves.md",
                "Utilities" => "utilities.md",
                "Extending" => "extending.md",
                "Similar tools" => "similartools.md"],
    sitename="PressureDrop.jl",
    format = Documenter.HTML(prettyurls = isCI)
)

#tag versions before pushing
if isCI
    deploydocs(repo = "github.com/jnoynaert/PressureDrop.jl.git")
end