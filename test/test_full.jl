#%% test segment calculation - B&B vertical

pressurecorrelation = BeggsAndBrill
p_initial = 346.6
dp_est = 100
t_avg = 107.2
md_initial = 0
md_end = 3700
tvd_initial = 0
tvd_end = 3700
inclination = 0
id = 2.441
roughness = 0.0003
q_o = 375
q_w = 0
GLR = 480
APIoil = 41.06
sg_water = 1
sg_gas = 0.916
molFracCO2 = 0
molFracH2S = 0
pseudocrit_pressure_correlation = HankinsonWithWichertPseudoCriticalPressure
pseudocrit_temp_correlation = HankinsonWithWichertPseudoCriticalTemp
Z_correlation = PapayZFactor
gas_viscosity_correlation = LeeGasViscosity
solutionGORcorrelation = StandingSolutionGOR
oilVolumeFactor_correlation = StandingOilVolumeFactor
waterVolumeFactor_correlation = GouldWaterVolumeFactor
dead_oil_viscosity_correlation = GlasoDeadOilViscosity
live_oil_viscosity_correlation = ChewAndConnallySaturatedOilViscosity
error_tolerance = 1.0 #psi

dpx = calculate_pressuresegment_topdown(pressurecorrelation, p_initial, dp_est, t_avg,
                                md_initial, md_end, tvd_initial, tvd_end, inclination, id, roughness,
                                q_o, q_w, GLR, APIoil, sg_water, sg_gas, molFracCO2, molFracH2S,
                                pseudocrit_pressure_correlation, pseudocrit_temp_correlation, Z_correlation,
                                gas_viscosity_correlation, solutionGORcorrelation, oilVolumeFactor_correlation, waterVolumeFactor_correlation,
                                dead_oil_viscosity_correlation, live_oil_viscosity_correlation, error_tolerance)

@test dpx ≈ 3770 * (0.170 + 0.312) / 2 atol = 100 #order of magnitude test based on average between example points in Takacs (50).


#TODO: extensive segment and integration testing for each correlation; minor errors may be propogating but not detected!
