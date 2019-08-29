# PressureDrop.jl Documentation

```@meta
CurrentModule = PressureDrop
```

PressureDrop.jl is a Julia package for computing multiphase pressure profiles for gas lift optimization of oil & gas wells.

Outlet-referenced models for producing wells using non-coupled temperature gradients are currently supported. 

!!! note

    Note that all inputs and calculations are in U.S. field units.

# Overview

The pressure traverse along the producing flow path of an oil well (or its bottomhole pressure at the interface between the well and the reservoir) is critical information for many production and reservoir engineering workflows. This can be expensive to measure directly in many conditions, so 1D pressure and temperature models are frequently used to infer producing conditions.

In most wells, the fluid flow contains three distinct phases (oil, water, and gas) with temperature- and pressure-dependent properties, including density, viscosity. In addition, the fluids will exhibit varying degrees of miscibility and entrainment depending on pressure and temperature as well, such that the pressure change with respect to distance or depth along the wellbore is itself dependent on pressure:

``\frac{∂p}{∂h} = f(p,T)``

Since the temperature profile varies according to depth, when assuming steady-state flow, consistent composition of each fluid, and a steady-state temperature profile uncoupled from pressure, the above can be re-expressed in terms of distance and pressure only:

``\frac{dp}{dh} = f(p,h)``

Note that ``f`` is composited from many empirical functions¹ and in most cases cannot be expressed in a tractable analytical form when dealing with multiphase flow. This is further complicated in gas lift wells, where the injection point is also pressure dependent, but changes the fluid composition and pressure profile above it.

Currently, no mechanistic or empirical correlations fully capture the variability in these three-phase fluid flows (or the direct properties of the fluids themselves), so the most performant methods are typically matched to the conditions and fluids they were developed to describe².

Most techniques for applying these correlations to calculate pressure profiles involve dividing the wellbore into a 1-dimensional series of discrete segments and calculating an average pressure change for the entire segment. Increasing segmentation will generally improve accuracy at the cost of additional computation. The most feasible and stable method for resolving the pressure change in each segment is typically some form of fixed-point iteration by specifying an error tolerance ε and iterating until ``f(p) - p < ε``, where p is the average pressure in the segment.

---

¹For an example of many widely-used correlations for pressure- and temperature-dependent properties of oil, water, and gas, see the "Fluids" section of the [Fekete documentation](http://www.fekete.com/SAN/TheoryAndEquations/HarmonyTheoryEquations/Content/HTML_Files/Reference_Material/Calculations_and_Correlations/Calculations_and_Correlations.htm).

²For a comprehensive overview of the theory of these applied methods in the context of gas lift engineering and nodal analysis, see Gábor Takács' *Gas Lift Manual* (2005, Pennwell Books).

# Contents

```@contents
Pages = [
  "core.md",
  "utilities.md",
  "plotting.md",
  "correlations.md",
  "pvt.md",
  "valves.md",
  "utilities.md",
  "extending.md"
]
Depth = 2
```
