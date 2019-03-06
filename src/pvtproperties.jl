# Helper functions for pvt properties.

#%% Gas

"""
Gas viscosity (μ_g) in centipoise.

Takes gas specific gravity, psia, °F, Z (deviation factor).

Lee et al. Takacs p20.
https://petrowiki.org/Gas_viscosity
"""
function LeeGasViscosity(specificGravity, psiAbs, tempF, Z)

  tempR = tempF + 459.67
  molecularWeight = 28.967 * specificGravity
  density = (psiAbs * molecularWeight * 0.00149406) / (Z*tempR)

  K = ((0.00094 + 2.0* 10^-6.0 *molecularWeight) * tempR^1.5 ) / (200.0 + 19.0*molecularWeight + tempR)
  X = 3.5 + (986.0/tempR) + 0.01*molecularWeight
  Y = 2.4 - 0.2*X

  return K * exp(X * density^Y ); #in centipoise
end


"""
Pseudo-critical temperature, adjusted pseudo-critical temperature in °R, and Wichert correction factor

Takes gas specific gravity, mol fraction of CO₂, mol fraction of H₂S.

Hankin-Thomas-Phillips, with Wichert and Aziz correction for sour components. Takacs p18.
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
Pseudo-critical pressure in psia.

Takes gas specific gravity, mol fraction of CO₂, mol fraction of H₂S.

Hankin-Thomas-Phillips, with Wichert and Aziz correction for sour components. Takacs p18.
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
Natural gas deviation factor (Z).

Take pseudocritical pressure (psia), pseudocritical temperature (°R), pressure (psia), temperature (°F).

Papay. Takacs p19.
"""
function PapayZFactor(pressurePseudoCritical, tempPseudoCriticalRankine, psiAbs, tempF)

  pressurePseudoCriticalReduced = psiAbs/pressurePseudoCritical

  tempPseudoCriticalReduced = (tempF + 459.67)/tempPseudoCriticalRankine

  return 1 - (3.52*pressurePseudoCriticalReduced)/(10^ (0.9813*tempPseudoCriticalReduced)) +
  (0.274*pressurePseudoCriticalReduced*pressurePseudoCriticalReduced)/(10^(0.8157*tempPseudoCriticalReduced))
end


"""
Corrected gas volume factor (B_g).

Takes absolute pressure (psia), Z-factor, temp in °F.
"""
function gasVolumeFactor(pressureAbs, Z, tempF)

  return 0.0283 * (Z*(tempF+459.67)/pressureAbs)
end



#%% Oil
#TODO: add Hanafy et al from http://fekete.com/SAN/TheoryAndEquations/WellTestTheoryEquations/Hanafy.htm
#TODO: add Vasquez-Beggs

"""
Solution GOR (Rₛ) in scf/bbl.

Takes oil gravity (°API), gas specific gravity, pressure (psia), temp (°F).

Standing. Takacs p13.
"""
function StandingSolutionGOR(APIoil, specificGravityGas, psiAbs, tempF)

  y = 0.00091*tempF - 0.0125*APIoil

  return specificGravityGas * (psiAbs / (18 * 10^y) )^1.205 #scf/bbl
end


"""
Standing oil volume factor (Bₒ).

Takes oil gravity (°API), gas specific gravity, solution GOR (scf/bbl), absolute pressure (psia), temp (°F).

Standing. Takacs p13.
"""
function StandingOilVolumeFactor(APIoil, specificGravityGas, solutionGOR, psiAbs, tempF)

  fFactor = solutionGOR * sqrt(specificGravityGas/(141.5/(APIoil + 131.5))) + 1.25*tempF

  return 0.972 + (0.000147) * fFactor^1.175
end


"""
Oil density (ρₒ) in mass-lbs per ft³.

Takes oil gravity (°API), gas specific gravity, solution GOR (scf/bbl), oil volume factor.
"""
function oilDensity_insitu(APIoil,  specificGravityGas,  solutionGOR,  oilVolumeFactor)

  return (141.5/(APIoil + 131.5)*62.42796 + 0.0136*specificGravityGas*solutionGOR)/oilVolumeFactor #mass-lbs per ft³
end


"""
Dead oil viscosity (μ_oD) in centipoise.

Takes oil gravity (°API), temp (°F).

Use with caution at 100-150° F: viscosity can be significantly overstated.

Beggs and Robinson.
"""
function BeggsAndRobinsonDeadOilViscosity(APIoil,  tempF)

  if 100 <= tempF <= 155
    print("Warning: using Beggs and Robinson for dead oil viscosity at $tempF--consider using another correlation.")
  end

  y = 10^(3.0324 - 0.02023*APIoil)
  x = y*tempF^-1.163

  return 10^x - 1
end


"""
Dead oil viscosity (μ_oD) in centipoise.

Takes oil gravity (°API), temp (°F).

Glaso. https://petrowiki.org/Calculating_PVT_properties#Dead_oil_viscosity
"""
function GlasoDeadOilViscosity(APIoil,  tempF)

  return ((3.141* 10^10) / tempF^3.444) * log10(APIoil)^(10.313*log10(tempF) - 36.447)
end


"""
Saturated oil viscosity (μₒ) in centipoise.

Takes dead oil viscosity (cp), solution GOR (scf/bbl).

Chew and Connally method. Takacs p15.
"""
function ChewAndConnallySaturatedOilViscosity(deadOilViscosity,  solutionGOR)

   A = 0.2 + 0.8 / 10.0^(0.00081*solutionGOR)
   b = 0.43 + 0.57 / 10.0^(0.00072*solutionGOR)

  return A * deadOilViscosity^b
end


#%% Water

"""
Water density in lb per ft³.

Takes gravity.
"""
function waterDensity_stb(waterGravity)
  return waterGravity * 62.4 #lb per ft^3
end


"""
Water volume factor (B_w).

Takes absolute pressure (psia), temp (°F).

Gould. Takacs p12.
"""
function GouldWaterVolumeFactor(pressureAbs,  tempF)

  return 1.0 + 1.21 * 10^-4 * (tempF-60)  + 10^-6 * (tempF-60)^2 - 3.33 * pressureAbs * 10^-6
end

const assumedWaterViscosity = 1.0 #centipoise
