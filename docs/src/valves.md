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
valves = GasliftValves([1813,2375,2885,3395], #valve MDs
                       [1005,990,975,960], #valve PTROs (psig)
                       [0.073,0.073,0.073,0.073], #valve R-values
                       [16,16,16,16]) #valve port sizes in 64ths inches

tubing_pressures = [150,837,850,840,831] #pressures at depth
casing_pressures = 1070 .+ [0,53,70,85,100]
temps = [135,145,148,151,153] #temps at depth

vdata, inj_depth = valve_calcs(valves = valves, well = well, sg_gas = 0.72, tubing_pressures = tubing_pressures, casing_pressures = casing_pressures, tubing_temps = temps, casing_temps = temps)

valve_table(vdata, inj_depth)
```

## Functions

```@autodocs
Modules = [PressureDrop]
Pages   = ["valvecalculations.jl"]
```
