# Helper functions for pvt properties.

#%% Gas

"""
`LeeGasViscosity(specificGravity, psiAbs, tempF, Z)`

Gas viscosity (μ_g) in centipoise.

Takes gas specific gravity, psia, °F, Z (deviation factor).

Lee et al 1966 method.
"""
function LeeGasViscosity(specificGravity, psiAbs, tempF, Z)

  tempR = tempF + 459.67
  molecularWeight = 28.967 * specificGravity
  density = (psiAbs * molecularWeight * 0.00149406) / (Z*tempR)

  K = ((0.00094 + 2.0* 10^-6.0 *molecularWeight) * tempR^1.5 ) / (200.0 + 19.0*molecularWeight + tempR)
  X = 3.5 + (986.0/tempR) + 0.01*molecularWeight
  Y = 2.4 - 0.2*X

  return K * exp(X * density^Y ) #in centipoise
end


"""
`HankinsonWithWichertPseudoCriticalTemp(specificGravity, molFracCO2, molFracH2S)`

Pseudo-critical temperature, adjusted pseudo-critical temperature in °R, and Wichert correction factor

Takes gas specific gravity, mol fraction of CO₂, mol fraction of H₂S.

Hankin-Thomas-Phillips method, with Wichert and Aziz correction for sour components.
"""
function HankinsonWithWichertPseudoCriticalTemp(specificGravity, molFracCO2, molFracH2S)
  #Hankinson-Thomas-Phillips pseudo-parameter:
  tempPseudoCritical = 170.5 + (307.3*specificGravity)

  #Wichert and Aziz pseudo-param correction for sour gas components:
  A = molFracCO2 + molFracH2S
  B = molFracH2S
  correctionFactor = 120.0*(A^0.9 - A^1.6) + 15.0*(B^0.5 - B^4.0)

  tempPseudoCriticalMod = tempPseudoCritical - correctionFactor

  return tempPseudoCritical, tempPseudoCriticalMod, correctionFactor
end


"""
`HankinsonWithWichertPseudoCriticalPressure(specificGravity, molFracCO2, molFracH2S)`

Pseudo-critical pressure in psia.

Takes gas specific gravity, mol fraction of CO₂, mol fraction of H₂S.

Hankin-Thomas-Phillips method, with Wichert and Aziz correction for sour components.
"""
function HankinsonWithWichertPseudoCriticalPressure(specificGravity, molFracCO2, molFracH2S)
  #Hankinson-Thomas-Phillips pseudo-parameter:
  pressurePseudoCritical = 709.6 - (58.7*specificGravity) #I think the 58.7 is incorrectly identified as 56.7 in Takacs text

  #pseudo-temperatures:
  tempPseudoCritical, tempPseudoCriticalMod, correctionFactor = HankinsonWithWichertPseudoCriticalTemp(specificGravity, molFracCO2, molFracH2S)
  B = molFracH2S

  return pressurePseudoCritical * tempPseudoCriticalMod / (tempPseudoCritical + B*(1-B)*correctionFactor)
end


"""
`PapayZFactor(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)`

Natural gas compressibility deviation factor (Z).

Take pseudocritical pressure (psia), pseudocritical temperature (°R), pressure (psia), temperature (°F).

Papay 1968 method.
"""
function PapayZFactor(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)

  pressurePseudoCriticalReduced = psiAbs/pressurePseudoCritical

  tempPseudoCriticalReduced = (tempF + 459.67)/tempPseudoCriticalRankine

  return 1 - (3.52*pressurePseudoCriticalReduced)/(10^ (0.9813*tempPseudoCriticalReduced)) +
  (0.274*pressurePseudoCriticalReduced*pressurePseudoCriticalReduced)/(10^(0.8157*tempPseudoCriticalReduced))
end


"""
`KareemEtAlZFactor(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)`

Natural gas compressibility deviation factor (Z).

Take pseudocritical pressure (psia), pseudocritical temperature (°R), pressure (psia), temperature (°F).

Can use gauge pressures so long as unit basis matches.

Direct correlation continuous over 0.2 ≤ P_pr ≤ 15.

Kareem, Iwalewa, Al-Marhoun, 2016.
"""
function KareemEtAlZFactor(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)

  P_pr = psiAbs/pressurePseudoCritical

  T_pr = (tempF + 459.67)/tempPseudoCriticalRankine

  if !(0.2 ≤ P_pr ≤ 15) | !(1.15 ≤ T_pr ≤ 3)
    @info "Using Kareem et al Z-factor correlation with values outside of 0.2 ≤ P_pr ≤ 15, 1.15 ≤ T_pr ≤ 3."
  end

  a1 = 0.317842
  a2 = 0.382216
  a3 = -7.768354
  a4 = 14.290531
  a5 = 0.000002
  a6 = -0.004693
  a7 = 0.096254
  a8 = 0.16672
  a9 = 0.96691
  a10 = 0.063069
  a11 = -1.966847
  a12 = 21.0581
  a13 = -27.0246
  a14 = 16.23
  a15 = 207.783
  a16 = -488.161
  a17 = 176.29
  a18 = 1.88453
  a19 = 3.05921

  t = 1 / T_pr
  A = a1 * t * exp(a2 * (1-t)^2) * P_pr
  B = a3 * t + a4 * t^2 + a5 * t^6 * P_pr^6
  C = a9 + a8 * t * P_pr + a7 * t^2 * P_pr^2 + a6 * t^3 * P_pr^3
  D = a10 * t * exp(a11 * (1-t)^2)
  E = a12 * t + a13 * t^2 + a14 * t^3
  F = a15 * t + a16 * t^2 + a17 * t^3
  G = a18 + a19 * t

  y = D * P_pr / ((1 + A^2) / C - (A^2 * B) / C^3)
  z = D * P_pr * (1 + y + y^2 - y^3) / (D * P_pr + E * y^2 - F * y^G) / (1 - y)^3

  return z
end


"""
`KareemEtAlZFactor_simplified(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)`

Natural gas compressibility deviation factor (Z).

Take pseudocritical pressure (psia), pseudocritical temperature (°R), pressure (psia), temperature (°F).

Linearized form from Kareem, Iwalewa, Al-Marhoun, 2016.
"""
function KareemEtAlZFactor_simplified(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)

  P_pr = psiAbs/pressurePseudoCritical

  T_pr = (tempF + 459.67)/tempPseudoCriticalRankine

  if !(0.2 ≤ P_pr ≤ 15) | !(1.15 ≤ T_pr ≤ 3)
    @info "Using Kareem et al Z-factor correlation with values outside of 0.2 ≤ P_pr ≤ 15, 1.15 ≤ T_pr ≤ 3."
  end

  a1 = 0.317842
  a2 = 0.382216
  a3 = -7.768354
  a4 = 14.290531
  a5 = 0.000002
  a6 = -0.004693
  a7 = 0.096254
  a8 = 0.16672
  a9 = 0.96691


  t = 1 / T_pr
  A = a1 * t * exp(a2 * (1-t)^2) * P_pr
  B = a3 * t + a4 * t^2 + a5 * t^6 * P_pr^6
  C = a9 + a8 * t * P_pr + a7 * t^2 * P_pr^2 + a6 * t^3 * P_pr^3

  z = (1+ A^2)/C - A^2 * B / C^3

  return z
end


"""
`gasVolumeFactor(pressureAbs, Z, tempF)`

Corrected gas volume factor (B_g).

Takes absolute pressure (psia), Z-factor, temp (°F).
"""
function gasVolumeFactor(pressureAbs, Z, tempF)

  return 0.0283 * (Z*(tempF+459.67)/pressureAbs)
end


"""
`gasDensity_insitu(specificGravityGas, Z_factor, abspressure, tempF)`

In-situ gas density in lb/ft³ (ρ_g).

Takes gas s.g., Z-factor, absolute pressure (psia), temperature (°F).
"""
function gasDensity_insitu(specificGravityGas, Z_factor, abspressure, tempF)

  tempR = tempF + 459.67

  return 2.7 * specificGravityGas * abspressure / (Z_factor * tempR)
end


#%% Oil
#TODO: add Vasquez-Beggs for solution GOR
#TODO: add Hanafy et al correlations

"""
`StandingBubblePoint(APIoil, sg_gas, R_b, tempF)`

Bubble point pressure in psia.

Takes oil gravity (°API), gas specific gravity, total solution GOR at pressures above bubble point (R_b, scf/bbl), temp (°F).
"""
function StandingBubblePoint(APIoil, sg_gas, R_b, tempF)
    y = 0.00091 * tempF - 0.0125 * APIoil
    return 18.2 * ((R_b/sg_gas)^0.83 * 10^y - 1.4)
end


"""
`StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF)`

Solution GOR (Rₛ) in scf/bbl.

Takes oil gravity (°API), gas specific gravity, pressure (psia), temp (°F), total solution GOR (R_b, scf/bbl), and bubblepoint value (psia).

Standing method.
"""
function StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF, R_b, bubblepoint::Real)

    if psiAbs >= bubblepoint
        return R_b
    else
        y = 0.00091*tempF - 0.0125*APIoil

        return specificGravityGas * (psiAbs / (18 * 10^y) )^1.205 #scf/bbl
    end
end


"""
`StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF)`

Solution GOR (Rₛ) in scf/bbl.

Takes oil gravity (°API), gas specific gravity, pressure (psia), temp (°F), total solution GOR (R_b, scf/bbl), and bubblepoint function.

Standing method.
"""
function StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF, R_b, bubblepoint::Function)

    if psiAbs >= bubblepoint(APIoil, specificGravityGas, R_b, tempF)
        return R_b
    else
        y = 0.00091*tempF - 0.0125*APIoil

        return specificGravityGas * (psiAbs / (18 * 10^y) )^1.205 #scf/bbl
    end
end


"""
`StandingOilVolumeFactor(APIoil, specificGravityGas, solutionGOR, psiAbs, tempF)`

Oil volume factor (Bₒ).

Takes oil gravity (°API), gas specific gravity, solution GOR (scf/bbl), absolute pressure (psia), temp (°F).

Standing method.
"""
function StandingOilVolumeFactor(APIoil, specificGravityGas, solutionGOR, psiAbs, tempF)

  fFactor = solutionGOR * sqrt(specificGravityGas/(141.5/(APIoil + 131.5))) + 1.25*tempF

  return 0.972 + (0.000147) * fFactor^1.175
end


"""
`oilDensity_insitu(APIoil,  specificGravityGas,  solutionGOR,  oilVolumeFactor)`

Oil density (ρₒ) in mass-lbs per ft³.

Takes oil gravity (°API), gas specific gravity, solution GOR (R_s, scf/bbl), oil volume factor.
"""
function oilDensity_insitu(APIoil,  specificGravityGas,  solutionGOR,  oilVolumeFactor)

  return (141.5/(APIoil + 131.5)*350.4 + 0.0764*specificGravityGas*solutionGOR)/(5.61 * oilVolumeFactor) #mass-lbs per ft³
end


"""
`BeggsAndRobinsonDeadOilViscosity(APIoil, tempF)`

Dead oil viscosity (μ_oD) in centipoise.

Takes oil gravity (°API), temp (°F).

Use with caution at 100-150° F: viscosity can be significantly overstated.

Beggs and Robinson method.
"""
function BeggsAndRobinsonDeadOilViscosity(APIoil, tempF)

  if 100 <= tempF <= 155
    @info "Warning: using Beggs and Robinson for dead oil viscosity at $(tempF)° F--consider using another correlation for 100-150° F."
  end

  y = 10^(3.0324 - 0.02023*APIoil)
  x = y*tempF^-1.163

  return 10^x - 1
end


"""
`GlasoDeadOilViscosity(APIoil, tempF)`

Dead oil viscosity (μ_oD) in centipoise.

Takes oil gravity (°API), temp (°F).

Glaso method.
"""
function GlasoDeadOilViscosity(APIoil, tempF)

  return ((3.141* 10^10) / tempF^3.444) * log10(APIoil)^(10.313*log10(tempF) - 36.447)
end


"""
`ChewAndConnallySaturatedOilViscosity(deadOilViscosity,  solutionGOR)`

Saturated oil viscosity (μₒ) in centipoise.

Takes dead oil viscosity (cp), solution GOR (scf/bbl).

Chew and Connally method to correct from dead to live oil viscosity.
"""
function ChewAndConnallySaturatedOilViscosity(deadOilViscosity,  solutionGOR)

   A = 0.2 + 0.8 / 10.0^(0.00081*solutionGOR)
   b = 0.43 + 0.57 / 10.0^(0.00072*solutionGOR)

  return A * deadOilViscosity^b
end


#%% Water

"""
`waterDensity_stb(waterGravity)`

Water density in lb per ft³.

Takes water specific gravity.
"""
function waterDensity_stb(waterGravity)
  return waterGravity * 62.4 #lb per ft^3
end


"""
`waterDensity_insitu(waterGravity, B_w)`

Water density in lb per ft³.

Takes specific gravity, B_w.
"""
function waterDensity_insitu(waterGravity, B_w)
  return waterGravity * 62.4 / B_w #lb per ft^3
end


"""
`GouldWaterVolumeFactor(pressureAbs,  tempF)`

Water volume factor (B_w).

Takes absolute pressure (psia), temp (°F).

Gould method.
"""
function GouldWaterVolumeFactor(pressureAbs,  tempF)

  return 1.0 + 1.21 * 10^-4 * (tempF-60)  + 10^-6 * (tempF-60)^2 - 3.33 * pressureAbs * 10^-6
end

const assumedWaterViscosity = 1.0 #centipoise


#%% Interfacial tension

"""
`gas_oil_interfacialtension(APIoil, pressureAbsolute, tempF)`

Gas-oil interfactial tension in dynes/cm.

Takes oil API, absolute pressure (psia), temp (°F).

Possibly Baker Swerdloff method; same method utilized by [Fekete](https://www.ihsenergy.ca/support/documentation_ca/Harmony/content/html_files/reference_material/calculations_and_correlations/pressure_loss_calculations.htm)
"""
function gas_oil_interfacialtension(APIoil, pressureAbsolute, tempF)

  if tempF <= 68.0
    σ_dead = 39 - 0.2571 * APIoil
  elseif tempF >= 100.0
    σ_dead = 37.5 - 0.2571 * APIoil
  else
    σ_68 = 39 - 0.2571 * APIoil
    σ_100 = 37.5 - 0.2571 * APIoil
    σ_dead = σ_68 + (tempF - 68.0) * (σ_100 - σ_68) / (100.0 - 68.0) # interpolate
  end

  C = 1.0 - 0.024 * pressureAbsolute^0.45
  return C * σ_dead
end


"""
`gas_water_interfacialtension(pressureAbsolute, tempF)`

Gas-water interfactial tension in dynes/cm.

Takes absolute pressure (psia), temp (°F).

Possibly Baker Swerdloff method; same method utilized by [Fekete](https://www.ihsenergy.ca/support/documentation_ca/Harmony/content/html_files/reference_material/calculations_and_correlations/pressure_loss_calculations.htm)
"""
function gas_water_interfacialtension(pressureAbsolute, tempF)

  if tempF <= 74
    return 75 - 1.018 * pressureAbsolute^0.349
  elseif tempF >= 280
    return 53 - 0.1048 * pressureAbsolute^0.637
  else
    σ_74 = 75 - 1.018 * pressureAbsolute^0.349
    σ_280 = 53 - 0.1048 * pressureAbsolute^0.637
    return σ_74 + (tempF - 74.0) * (σ_280 - σ_74) / (280.0 - 74.0) #interpolate
  end
end
