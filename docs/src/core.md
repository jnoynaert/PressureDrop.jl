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
                  temperature_method = "linear", #temperatures can be calculated or provided directly as a array
                  WHT = 105, BHT = 160, #°F
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
tubing_pressures = pressure_and_temp!(model) #note that this updates temperature in the .temperatureprofile field of the WellModel
```

Several [plotting functions](@ref Plotting) are available to visualize the outputs.

```@example core
using Gadfly #necessary to load plotting functions

plot_pressure(model, tubing_pressures, "Tubing Pressure Drop")
draw(SVG("plot-pressure-core.svg", 6inch, 4inch), ans); nothing # hide
```

![](plot-pressure-core.svg)

Pressure traverses for just tubing or just casing, utilizing an existing temperature profile, can be calculated using [`traverse_topdown`](@ref) or [`casing_traverse_topdown`](@ref).

### Gas lift analysis

The [`gaslift_model!`](@ref) function will calculate the pressure and temperature profiles, most likely operating point (assuming single-point injection), and opening and closing pressures of the valves.

```@example core
tubing_pressures, casing_pressures, valvedata = gaslift_model!(model, find_injectionpoint = true,
               dp_min = 100) #required minimum ΔP at depth to consider as an operating valve

plot_gaslift(model, tubing_pressures, casing_pressures, valvedata, "Gas Lift Analysis Plot")
draw(SVG("plot-gl-core.svg", 6inch, 4inch), ans); nothing # hide
```

![](plot-gl-core.svg)

The results of the valve calculations can be printed as a table:

```@example core
valve_table(valvedata)
```

The data for a valve table can be calculated directly using [`valve_calcs`](@ref), which will interpolate pressures and temperatures at depth from known producing P/T profiles.

## Types

```@docs
Wellbore
GasliftValves
WellModel
```

## Functions

```@docs
traverse_topdown
casing_traverse_topdown
pressure_and_temp!
pressures_and_temp!
gaslift_model!
```
