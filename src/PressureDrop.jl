module PressureDrop

# Calculate pressure and temperature profiles for oil & gas wells.

#TODO: avoid duplicating calculations for no reason
#TODO: Finish profiling using @time, @code_warntype, Traceur, ProfileView
#TODO: @inbounds master loop
#TODO: add test/runtests.jl per Takacs results, or per a reference result from IHS, and run on every new build once you have a main loop in place

#using
#import


include("pvtproperties.jl")
include("pressurecorrelations.jl")
include("tempcorrelations.jl")



"""
Develop pressure traverse from wellhead down to datum.

Inputs: survey, WHP, WHT, BHT,
GOR,
gas s.g., water s.g., oil API
which correlations to use for fluid properties
which pressure drop correlation to use

Returns: pressure profile, temperature profile as two separate Array{Float64, 1}s.

Automatically meshes the deviation profile without resizing etc -- keep that separate as a utility.
"""
function traverse_topdown(;pressurecorrelation = BeggsAndBrill)
    # Takacs 113 figure 2-41
    # Make sure this is from wellhead down
    # TODO: calculate BHT?
    # TODO: use GOR to initialize model, since GLR is ambiguous term for gas lift opt, but GOR is understand to be reservoir #s
    # TODO: make sure to utilize named arguments

    set ϵ or error maximum
    set h_initial
    set dh
    estimate ΔP_est

    iteration 0:
    p_avg = p_initial + dp_est/2

    t_avg = t(h_initial + dh/2)

    calculate pvt properties at p_avg, t_avg

    calculate dp/dh #use pressurecorrelation() directly instead of pattern matching

    dp_calc = dh*(dp/dh)

    if abs(dp_est - dp_calc) > ϵ:
        new dp_est = dp_calc
        reiterate
    else
        h_initial_new = h_initial + dh
        p_initial_new = p_initial + dp_calc

    if h_initial_new >= td_well:
        proportionally estimate pwf based on dp_calc/dh
end

end #module PressureDrop
