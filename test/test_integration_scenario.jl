#General parameters
BHT = 165

WHP = 220 PSIG

oil_API = 35
sg_gas = 0.65
CO2 = 0.5%
N2 = 2%
sg_water = 1.07

# settings for PVT in IHS:
dead oil visco -> Glaso
Saturated oil vis -> chew & ChewAndConnallySaturatedOilViscosity
undersaturated -> vazquez
gas -> Lee
Water -> matthews & russell
oil density -> standing
bpp & R_S -> standing
oil comp -> vazquez & beggs
oil fvf -> standing
Z -> Hall & Yarborough

survey = Sawgrass 9 to TD (12175 MD), first segment edited to 0 inclination (460 md/tvd)

tubing id = 2.441

#%% Scenario A
500 bpd
WC = 50%
GLR = 4500
WHT = 100

#%% Scenario B
250 bpd
WC = 25%
GLR = 6000
WHT = 90

#%% Scenario C
1000 bpd
WC = 75%
GLR = 3000
WHT = 105

#%% Scenario D
3000 bpd
WC = 85%
GLR = 1200
WHT = 115

#%% Scenario E
50 bpd
WC = 25%
GLR = 10000
WHT = 80




# repeat for H&B mod--with and without Griffith and Wallis, B&B mod, Duns & Ros, Gray, and a mechanistic correlation

# check across the whole wellbore by interpolating
