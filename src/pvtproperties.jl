# Helper functions for pvt properties.

"""
Gas viscosity in centipoise.

Takes gas specific gravity, psia, °F, Z (deviation factor).

Lee et al. Takacs p20.
https://petrowiki.org/Gas_viscosity
"""
function LeeGasViscosity(specificGravity, psiAbs, tempF, Z)

  tempR = tempF + 459.67
  molecularWeight = 28.967 * specificGravity
  density = (psiAbs * molecularWeight * 0.00149406) / (Z*tempR)

  K1 = ( (0.00094 + 2.0* 10^-6.0 *molecularWeight) * tempR^1.5 ) / (200.0 + 19.0*molecularWeight + tempR)
  X = 3.5 + (986.0/tempR) + 0.01*molecularWeight
  Y = 2.4 - 0.2*X

  return K1 * exp( X * density^Y ); #in centipoise
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


## <-- resume work here
"""
"""
double PapayZFactor (double pressurePseudoCritical, double tempPseudoCriticalRankine, double psiAbs, double tempF) {

  double pressurePseudoCriticalReduced = psiAbs/pressurePseudoCritical;

  double tempPseudoCriticalReduced = (tempF + 459.67)/tempPseudoCriticalRankine;

  #Papay Z-factor correlation. Useful because it does not require solving iteratively
  return 1 - (3.52*pressurePseudoCriticalReduced)/(std::pow(10, 0.9813*tempPseudoCriticalReduced)) +
    (0.274*pressurePseudoCriticalReduced*pressurePseudoCriticalReduced)/std::pow(10, 0.8157*tempPseudoCriticalReduced);
}

# [[Rcpp::export]]
double GasVolumeFactor (double pressureAbs, double Z, double tempF) {
  return 0.0283 * (Z*(tempF+459.67)/pressureAbs);
}

#  Rcout << "pressurePseudoCritical: " << pressurePseudoCritical << "\n" ;


#**Oil PVT functions
#alternative: Hanafy et al http:#fekete.com/SAN/TheoryAndEquations/WellTestTheoryEquations/Hanafy.htm


double StandingSolutionGOR (double APIoil, double specificGravityGas, double psiAbs, double tempF) {  #only use

  double y = 0.00091*tempF - 0.0125*APIoil;
  return specificGravityGas * std::pow( psiAbs / (18 * std::pow(10, y)), 1.205); #results in scf/bbl
}

#volume factor

double StandingOilVolumeFactor (double APIoil, double specificGravityGas, double solutionGOR, double psiAbs, double tempF) {

  double fFactor = solutionGOR * std::sqrt(specificGravityGas/(141.5/(APIoil + 131.5))) + 1.25*tempF;

  return 0.972 + (0.000147)*std::pow(fFactor, 1.175);
}

#Oil density
double oilDensity (double APIoil, double specificGravityGas, double solutionGOR, double oilVolumeFactor) {

  return (141.5/(APIoil + 131.5)*62.42796 + 0.0136*specificGravityGas*solutionGOR)/oilVolumeFactor; #mass-pounds over cubic feet
}

#dead oil viscosity
double BeggsAndRobinsonDeadOilViscosity (double APIoil, double tempF) { #careful! this overstates viscosity at 100 to 150 F

  double specificGravity = 141.5/(APIoil + 131.5);
  double x = std::pow(tempF, -1.163) * std::exp(13.108 - (6.591/specificGravity));

  return std::pow(10, x) - 1;
}

double GlasoDeadOilViscosity (double APIoil, double tempF) {

  return ((3.141* std::pow(10.0,10.0)) / std::pow(tempF, 3.444)) * std::pow( std::log10(APIoil), 10.313*std::log10(tempF) - 36.447);

}


double ChewAndConnallySaturatedOilViscosity (double deadOilViscosity, double solutionGOR) {

  double A = 0.2 + 0.8/std::pow(10.0, 0.00081*solutionGOR);
  double B = 0.43 + 0.57/std::pow(10.0, 0.00072*solutionGOR);

  return A*std::pow(deadOilViscosity, B);
}


#**Water PVT functions
double waterDensity (double waterGravity) {
  return waterGravity * 62.4; #lb per ft^3
}


double GouldWaterVolumeFactor (double pressureAbs, double tempF) {

  return 1.0 + 1.21*std::pow(10, -4)*(tempF-60) + std::pow(10, -6)*std::pow(tempF-60, 2) - 3.33*pressureAbs*std::pow(10,-6);
}

CONST assumedWaterViscosity = 1.0 #centipoise
