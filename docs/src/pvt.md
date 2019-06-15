# PVT properties

Most of the PVT property functions are passed as arguments to the core calculation functions, but are not used directly.

### Exported functions

User-facing functions used for model construction.

#### Oil

- [`StandingSolutionGOR`](@ref)
- [`BeggsAndRobinsonDeadOilViscosity`](@ref)
- [`GlasoDeadOilViscosity`](@ref)
- [`ChewAndConnallySaturatedOilViscosity`](@ref)

#### Gas

- [`LeeGasViscosity`](@ref)
- [`HankinsonWithWichertPseudoCriticalTemp`](@ref)
- [`HankinsonWithWichertPseudoCriticalPressure`](@ref)
- [`PapayZFactor`](@ref)
- [`KareemEtAlZFactor`](@ref)
- [`KareemEtAlZFactor_simplified`](@ref)

#### Water

- [`GouldWaterVolumeFactor`](@ref)

### Internal functions

- [`PressureDrop.gasVolumeFactor`](@ref)
- [`PressureDrop.gasDensity_insitu`](@ref)
- [`PressureDrop.oilDensity_insitu`](@ref)
- [`PressureDrop.waterDensity_stb`](@ref)
- [`PressureDrop.waterDensity_insitu`](@ref)
- [`PressureDrop.gas_oil_interfacialtension`](@ref)
- [`PressureDrop.gas_water_interfacialtension`](@ref)

## Functions

```@autodocs
Modules = [PressureDrop]
Pages   = ["pvtproperties.jl"]
```

```@docs
PressureDrop.gasVolumeFactor
PressureDrop.gasDensity_insitu
PressureDrop.oilDensity_insitu
PressureDrop.waterDensity_stb
PressureDrop.waterDensity_insitu
PressureDrop.gas_oil_interfacialtension
PressureDrop.gas_water_interfacialtension
```
