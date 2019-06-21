# Correlations

Pressure, temperature, and friction factor correlations. Not used directly but passed as model arguments.

- Pressure drop correlations
    - [`Beggs and Brill`](@ref BeggsAndBrill), with Payne correction
    - [`Hagedorn and Brown`](@ref HagedornAndBrown), with Griffith and Wallis bubble flow correction

- Temperature correlations & methods
    - [`Linear temperature profile`](@ref linear_wellboretemp)
    - [`Shiu temperature profile`](@ref Shiu_wellboretemp): Ramey temperature correlation with Shiu relaxation factor
    - [`Ramey temp`](@ref Ramey_temp): single-point Ramey temperature correlation
    - [`Shiu relaxation factor`](@ref Shiu_Beggs_relaxationfactor)

- Friction factor correlations
    - [`Serghide friction factor`](@ref SerghideFrictionFactor)(preferred)
    - [`Chen friction factor`](@ref ChenFrictionFactor)

## Pressure correlations

```@docs
BeggsAndBrill
HagedornAndBrown
```

## Temperature correlations

```@docs
linear_wellboretemp
Shiu_wellboretemp
Ramey_temp
Shiu_Beggs_relaxationfactor
```

## Friction factor correlations

```@docs
SerghideFrictionFactor
ChenFrictionFactor
```
