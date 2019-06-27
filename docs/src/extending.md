# Extending the calculation engine

PVT or pressure/temperature correlations can easily be added, either by modifying the original source, or more simply by defining new functions in your script or session that match the interface of existing functions.

For example, to add a new PVT function, first inspect either the function signature, source, or documentation for one of the functions of the category you are adding to:
```
julia> using PressureDrop

help?> StandingSolutionGOR
```
```
StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF, R_b, bubblepoint::Real)

Solution GOR (Rₛ) in scf/bbl.

Takes oil gravity (°API), gas specific gravity, pressure (psia), temp (°F), total solution GOR (R_b, scf/bbl), and bubblepoint value (psia).

<other methods not shown>
```

Then simply define your new function, making sure to either utilize the same interface, or capture extra unneeded arguments.
```
function HanafySolutionGOR(APIoil, specificGravityGas, psiAbs, tempF, R_b, bubblepoint::Real)
  if

  elseif psiAbs <= 157.28
    return 0
  else
    return -49.069 + 0.312 * psiAbs
end
```

The new PVT function can now be added to an existing model (or used in a pressure traverse function call or the creation of a new model):
```
oldmodel.solutionGORcorrelation = HanafySolutionGOR
```

Note that in this example the new method for solution GOR will only handle bubble points defined in absolute pressure.

Utilizing `RCall`, `PyCall`, or `ccall` will also allow adding functions defined in R, Python, and C or Fortran respectively.

 As of v1.0, defining new correlations or model functions that require additional arguments is not supported without modifying the source. However, feel free to either submit a pull request or open an issue requesting the additional functionality (include reference to the source material and several test cases).
