# Plotting

All plotting functions wrap Gadfly plot definitions.

!!! note

    Plotting functionality is lazily loaded and not available until `Gadfly` has been loaded.

## Examples

```@setup plots
# outputs & inputs hidden when doc is generated

using PressureDrop

segments = 100
MDs = range(0, stop = 5000, length = segments) |> collect
incs = repeat([0], inner = segments)
TVDs = range(0, stop = 5000, length = segments) |> collect

well = Wellbore(MDs, incs, TVDs, 2.441)

testpath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata")
valvepath = joinpath(testpath, "valvedata_wrappers_1.csv")
valves = read_valves(path = valvepath, delim = ',', skiplines = 1) #implicit read_valves test

model = WellModel(wellbore = well, roughness = 0.0, valves = valves,
                    temperature_method = "Shiu", geothermal_gradient = 1.0, BHT = 200,
                    pressurecorrelation = HagedornAndBrown, WHP = 350 - pressure_atmospheric, dp_est = 25,
                    q_o = 100, q_w = 500, GLR = 1200, APIoil = 35, sg_water = 1.0, sg_gas = 0.8, CHP = 1000, naturalGLR = 0)

tubing_pressures, casing_pressures, valvedata = gaslift_model!(model, find_injectionpoint = true, dp_min = 100)
```

### [`plot_gaslift`](@ref)

```@example plots
using Gadfly

plot_gaslift(model, tubing_pressures, casing_pressures, valvedata, "Gas Lift Analysis Plot")
draw(SVG("plot-gl.svg", 5inch, 4inch), ans); nothing # hide
```

![](plot-gl.svg)

### [`plot_pressure`](@ref)

```@example plots

plot_pressure(model, tubing_pressures, "Tubing Pressure Drop")
draw(SVG("plot-pressure.svg", 4inch, 4inch), ans); nothing # hide
```

![](plot-pressure.svg)


### [`plot_pressures`](@ref)

```@example plots

plot_pressures(model, tubing_pressures, casing_pressures, "Tubing and Casing Pressures")
draw(SVG("plot-pressures.svg", 4inch, 4inch), ans); nothing # hide
```

![](plot-pressures.svg)


### [`plot_temperature`](@ref)

```@example plots

plot_temperature(model.wellbore, model.temperatureprofile, "Temperature Profile")
draw(SVG("plot-temperature.svg", 4inch, 4inch), ans); nothing # hide
```

![](plot-temperature.svg)


### [`plot_pressureandtemp`](@ref)

```@example plots

plot_pressureandtemp(model, tubing_pressures, casing_pressures, "Pressures and Temps")
draw(SVG("plot-pressureandtemp.svg", 5inch, 4inch), ans); nothing # hide
```

![](plot-pressureandtemp.svg)


## Functions

```@docs
plot_pressure
plot_pressures
plot_temperature
plot_pressureandtemp
plot_gaslift
```
