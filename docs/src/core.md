# Core functionality

```@contents
Pages = ["core.md"]
Depth = 3
```

```@setup core
using PressureDrop

surveyfilepath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata/Sawgrass_9_32/Test_survey_Sawgrass_9.csv")

valvefilepath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata/valvedata_wrappers_1.csv")
```

## Creating and updating models

Model definitions are created and stored as [`WellModel`](@ref) objects. Although all of the functionality of this package is exposed as pure functions, mutating and copying `WellModel`s is a much easier way to track and iterate on parameter sets.

### Wellbores

The key component required for the pressure drop calculations is a [`Wellbore`](@ref) object that defines the flow path in terms of directional survey points (measured depth, inclination, and true vertical depth) and tubular inner diameter.

`Wellbore` objects can be constructed from arrays, or from CSV files with `read_survey`, which includes some optional convenience arguments to change delimiters, skip lines, or truncate the survey. Tubing IDs do not have to be uniform and can be specified segment to segment.

```@example core
examplewell = read_survey(path = surveyfilepath, id = 2.441, maxdepth = 6500) #an outlet point at 0 MD is added if not present
```

### Valve designs

`GasliftValves` objects define the valve strings in terms of measured run depth, test rack opening pressure, R value (ratio of the area of the port to the area of the bellows), and port size.

These can also be constructed directly or from CSV files.

```@example core
examplevalves = read_valves(path = valvefilepath)
```

### Models & parameter sets

`WellModel`s do not have to be completely specified, but require defining the minimum fields for a simple pressure drop. In general, sensible defaults are selected for PVT functions. See the [documentation](@ref WellModel) for a list of optional fields.

Note that defining a valve string is optional if all that is desired is a normal pressure drop or temperature calculation.

```@example core
model = WellModel(wellbore = examplewell, roughness = 0.00065,
                  valves = examplevalves,
                  pressurecorrelation = BeggsAndBrill,
                  WHP = 200, #wellhead pressure, psig
                  CHP = 1050, #casing pressure, psig
                  dp_est = 25, #estimated ΔP by segment. Not critical
                  temperature_method = "Shiu", #temperatures can be calculated or provided directly as a array
                  BHT = 160, geothermal_gradient = 0.9,  #°F, °F/100'
                  q_o = 100, q_w = 500, #bpd
                  GLR = 2500, naturalGLR = 400, #scf/bbl
                  APIoil = 35, sg_water = 1.05, sg_gas = 0.65);
```
Printing a WellModel will display all of its defined and undefined fields.

An import aspect of model definitions is that they include the temperature profile. Passing a model object to a wrapper function that calculates both pressure and temperature will mutate the temperature profile associate with the model.

## Pressure & temperature calculations

### Pressure traverses & temperature profiles

Pressure and temperature profiles can be generated from a `WellModel` using [`pressure_and_temp!`](@ref) (for tubing calculations only) or [`pressures_and_temp!`](@ref) (to include casing calculations).

```@example core
tubing_pressures = pressure_and_temp!(model); #note that this updates temperature in the .temperatureprofile field of the WellModel
```

Several [plotting functions](@ref Plotting) are available to visualize the outputs.

```@example core
using Gadfly #necessary to load plotting functions

plot_pressure(model, tubing_pressures, "Tubing Pressure Drop")
draw(SVG("plot-pressure-core.svg", 4inch, 4inch), ans); nothing # hide
```

![](plot-pressure-core.svg)

Pressure traverses for just tubing or just casing, utilizing an existing temperature profile, can be calculated using [`traverse_topdown`](@ref) or [`casing_traverse_topdown`](@ref).

### Gas lift analysis

The [`gaslift_model!`](@ref) function will calculate the pressure and temperature profiles, most likely operating point (assuming single-point injection), and opening and closing pressures of the valves.

```@example core
tubing_pressures, casing_pressures, valvedata = gaslift_model!(model, find_injectionpoint = true,
               dp_min = 100) #required minimum ΔP at depth to consider as an operating valve

plot_gaslift(model, tubing_pressures, casing_pressures, valvedata, "Gas Lift Analysis Plot")
draw(SVG("plot-gl-core.svg", 5inch, 4inch), ans); nothing # hide
```

![](plot-gl-core.svg)

The results of the valve calculations can be printed as a table:

```@example core
valve_table(valvedata)
```

The data for a valve table can be calculated directly using [`valve_calcs`](@ref), which will interpolate pressures and temperatures at depth from known producing P/T profiles.

## Bulk calculations

Pressure drops can be calculated in bulk, either by passing model arguments to functions directly, or by mutating or copying model objects.

```@example core
nominal_rate(D_sei, b) = ((1-D_sei)^(-b) - 1)/b #secant decline rates to nominal rates, b ≠ 0
hyperbolic_rate(q_i, b, D_sei, t) = q_i / (1 + b * nominal_rate(D_sei, b) * t)^(1/b) #spot rate from a hyperbolic decline for t in years

# generate test data
q_i = 3000
b = 1.2
decline = 0.85
timesteps = range(0,2, step = 1/365)
declinedata = [hyperbolic_rate(q_i, b, decline, time) for time in timesteps]
noise = [randn() .* 15 for sample in timesteps]
testdata = max.(declinedata .+ noise, 0)

# check results
days = timesteps .* 365
plot(x = days, y = testdata, Geom.path,
     Guide.xlabel("Time (days)"),
     Guide.ylabel("Total Fluid (bpd)"),
     Scale.y_continuous(format = :plain, minvalue = 0))
draw(SVG("test-data.svg", 6inch, 4inch), ans); nothing # hide
```

![](test-data.svg)

```@example core
# set up and calculate pressure data
examplewell = read_survey(path = surveyfilepath, id = 2.441, maxdepth = 6500)

function timestep_pressure(rate, temp, watercut, GLR)
    temps = linear_wellboretemp(WHT = temp, BHT = 165, wellbore = examplewell)

    return traverse_topdown(wellbore = examplewell, roughness = 0.0065, temperatureprofile = temps,
                     pressurecorrelation = BeggsAndBrill, dp_est = 25, error_tolerance = 0.1,
                     q_o = rate * (1 - watercut), q_w = rate * watercut, GLR = GLR,
                     APIoil = 36, sg_water = 1.05, sg_gas = 0.65,
                     WHP = 120)[end]
end

wellhead_temps = range(125, 85, length = 731)
watercuts = range(1, 0.5, length = 731)
GLR = range(0, 5000, length = 731)

pressures = timestep_pressure.(testdata, wellhead_temps, watercuts, GLR)

# examine outputs
plot(x = days, y = pressures, Geom.path, Theme(default_color = "purple"),
     Guide.xlabel("Time (days)"),
     Guide.ylabel("Flowing Pressure (psig)"),
     Scale.y_continuous(format = :plain, minvalue = 0),
     Guide.title("FBHP Over Time"))
draw(SVG("pressure-data.svg", 6inch, 4inch), ans); nothing # hide
```

![](pressure-data.svg)

## Types and Functions

- Types
    - [`Wellbore`](@ref)
    - [`GasliftValves`](@ref)
    - [`WellModel`](@ref)
- Functions
    [`traverse_topdown`](@ref)
    [`casing_traverse_topdown`](@ref)
    [`pressure_and_temp!`](@ref)
    [`pressures_and_temp!`](@ref)
    [`gaslift_model!`](@ref)

### Types

```@docs
Wellbore
GasliftValves
WellModel
```

### Functions

```@docs
traverse_topdown
casing_traverse_topdown
pressure_and_temp!
pressures_and_temp!
gaslift_model!
```
