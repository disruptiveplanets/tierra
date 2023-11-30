from tierra import Target
from tierra.transmission import TransmissionSpectroscopy
from bisect import bisect
import matplotlib.pyplot as plt
import numpy as np
import os

#Read the data from the file
LCFile = "data/WASP39b_NIRSPEC_G395H.txt"
Wavelength, WavelengthBin, Depth, DepthErr = np.loadtxt(LCFile, unpack=True)
Depth*=1e6
DepthErr*=1e6
    
WavelengthLower = Wavelength-WavelengthBin
WavelengthUpper = Wavelength+WavelengthBin

#########################################################################
#########################################################################
#   plt.figure(figsize=(12,8))
#   plt.errorbar(Wavelength, Depth, yerr=DepthErr, capsize=3, marker="o")
#   plt.xlabel("Wavelength [microns]")
#   plt.ylabel("(Rp/Rs)$^2$")
#   plt.tight_layout()
#   plt.show()
#########################################################################
#########################################################################

StellarDict = {}
PlanetaryDict = {}


StellarDict['Mass']=0.93
StellarDict['MassErr']=0.09
StellarDict['Radius']=0.895
StellarDict['RadiusErr']=0.023



PlanetaryDict['Mass']=0.28*317.907
PlanetaryDict['MassErr']=0.03*317.907
PlanetaryDict['Radius']=1.27*11.2089
PlanetaryDict['RadiusErr']=0.04*11.2089

PlanetaryDict['PT']=1
PlanetaryDict['P0']=1.0         #This is the reference pressure at the 
PlanetaryDict['T0']=1252.0      #This is in atmosphere

PlanetaryDict['MR_CO2'] = 1e-4
PlanetaryDict['MR_CO'] = 1e-5
PlanetaryDict['MR_H2O'] = 1e-3
PlanetaryDict['MR_SO2'] = 1e-5
MR_H2 = (1-1e-4-1e-3-1e-5-1e-5)/1.15

print("The mixing ratio of hydrogen is given by:", MR_H2)
#Calculate the mixing ratio for hydrogen
PlanetaryDict['MR_H2'] = MR_H2
Molecule2Look = ["CO2", "CO", "H2O", "SO2", "H2"]


CurrentLocation = os.getcwd()
if "gridsan" in CurrentLocation:
    print("Running in MIT Supercloud")
    CS_BaseLocation =  "/home/gridsan/pniraula/CrossSection_JWST/CrossSectionTIERRA"
else:
    print("Running in the Desktop Aekalavya")
    CS_BaseLocation =  "/media/prajwal/LaCie1/CrossSectionFromSuperCloud4JWST/CrossSectionTIERRA"

#All the retrievals will be done with respect to HITRAN
CurrentSystem = Target.System(Molecule2Look, PlanetaryDict, StellarDict)
CurrentSystem.LoadCrossSection(CS_BaseLocation, SubFolder="CS_1", CIA_Flag=True)

#Check the impact of opacity 
#remove this later to make things more automatic
CurrentSystem.InitiateSystem(PlanetaryDict)
CurrentSystem.PT_Profile(zStep=0.25, ShowPlot=False)
T1 = TransmissionSpectroscopy(CurrentSystem, CIA=True)
T1.CalculateTransmission(CurrentSystem)
CurrentWavelength = CurrentSystem.WavelengthArray*1e4
CurrentModel = T1.Spectrum


colorName = ["red", "blue", "green", "orange", "brown", "cyan"]
#VelocityValues = [0, 10, 100, 1000]
VelocityValues = [0, 20, 40, 60, 80, 100]


plt.figure(figsize=(12,7))
for VCounter, V0 in enumerate(VelocityValues):
   ShiftedWavelength = (1.0+V0/299792.458)*CurrentWavelength
   BinnedModel = np.zeros(len(Wavelength))
   counter = 0
   for Wl, Wp in zip(WavelengthLower, WavelengthUpper):
       StartIndex = bisect(ShiftedWavelength, Wl)
       StopIndex = bisect(ShiftedWavelength, Wp)
       BinnedModel[counter] = np.mean(CurrentModel[StartIndex:StopIndex])
       counter+=1  

   if VCounter == 0:
       RefModel = BinnedModel 

   plt.subplot(211) 
   plt.plot(Wavelength, BinnedModel, label=str(V0), color=colorName[VCounter])

   if VCounter>0: 
      plt.subplot(212)
      plt.plot(Wavelength, (BinnedModel-RefModel)/DepthErr, label=str(V0), color=colorName[VCounter]) 
plt.subplot(211) 
plt.xlim(min(Wavelength), max(Wavelength))
plt.xlabel("Wavelength")
plt.ylabel("[Rp/Rs]^2")
plt.legend(loc=1)
plt.subplot(212) 
plt.xlabel("Wavelength")
plt.ylabel("Deviation [Sigma]")
plt.xlim(min(Wavelength), max(Wavelength))
plt.legend(loc=1)
plt.tight_layout()
plt.savefig("DopplerShift.png")
plt.savefig("DopplerShift.pdf")
plt.show()


