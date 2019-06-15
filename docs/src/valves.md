# Valve calculations

Functions to generate valve performance curves and valve pressure tables using current conditions.

```@index
Modules = [PressureDrop]
Pages = ["valves.md"]
```

## Example

```@example valves
using PressureDrop

MDs = [0,1813, 2375, 2885, 3395]
TVDs = [0,1800, 2350, 2850, 3350]
incs = [0,0,0,0,0]
id = 2.441

well = Wellbore(MDs, incs, TVDs, id)
valves = GasliftValves([1813,2375,2885,3395], [1005,990,975,960], [0.073,0.073,0.073,0.073], [16,16,16,16])

tubing_pressures = [150,837,850,840,831]
casing_pressures = 1070 .+ [0,53,70,85,100]
temps = [135,145,148,151,153]

vdata, inj_depth = valve_calcs(valves = valves, well = well, sg_gas = 0.72, tubing_pressures = tubing_pressures, casing_pressures = casing_pressures, tubing_temps = temps, casing_temps = temps)

valve_table(vdata, inj_depth)
```

## Functions

```@autodocs
Modules = [PressureDrop]
Pages   = ["valvecalculations.jl"]
```
