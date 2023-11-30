import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_left


def planck(L,T):
    '''
    L: Wavelength in microns
    T: Temperature in Kelvin
    Returns:
      Array of spectral density in a given wavelength and a temperature
    '''
    G = 6.672e-11 # m^3 kg^-1 s^-2, courtesy Wikipedia
    c = 299792458 # m s^-1, courtesy Wikipedia
    h = 6.626e-34 # J s, courtesy Wikipedia
    boltzmann = 1.38e-23 # J K^-1, courtesy Wikipedia
    L = L*1e-6 #convert microns to meters
    return 2.*h*(c**2.)/(L**5.)/np.expm1(h*c/(L*boltzmann*T))

Data = np.loadtxt("WASP_96.csv", delimiter=",")

Wavelength = Data[:,0]
WavelengthStep = np.diff(Wavelength)
WavelengthStep = np.concatenate(([WavelengthStep[0]],WavelengthStep))
WavelengthLower = Wavelength-WavelengthStep/2
WavelengthUpper = Wavelength+WavelengthStep/2


TransitDepth = Data[:,1]


#Scale the error 
ErrorValueRef = (14544-14410)/2
WavelengthRef = 1.3708
ClosestIndex = bisect_left(Wavelength, WavelengthRef)

print("The closest index is given by", ClosestIndex)
print("The error value is given by", ClosestIndex)

WavelengthHR = np.linspace(min(WavelengthLower)/2.0, max(WavelengthUpper)+5.0, 100000)
WASP_BlackBody = planck(WavelengthHR, 5540)


WavelengthLowerRef = WavelengthLower[ClosestIndex] 
WavelengthUpperRef = WavelengthUpper[ClosestIndex]
StartIndex = bisect_left(WavelengthHR, WavelengthLowerRef)
StopIndex = bisect_left(WavelengthHR, WavelengthUpperRef)
Ref_Photon = np.sum(WASP_BlackBody[StartIndex:StopIndex])

PhotonsCount = []
PhotonsMean = []
for WStart, WStop in zip(WavelengthLower, WavelengthUpper):
  StartIndex = bisect_left(WavelengthHR, WStart)
  StopIndex = bisect_left(WavelengthHR, WStop)
  PhotonsCount.append(np.sum(WASP_BlackBody[StartIndex:StopIndex]))
  PhotonsMean.append(np.mean(WASP_BlackBody[StartIndex:StopIndex]))
PhotonsCount = np.array(PhotonsCount)


TransitError = ErrorValueRef*np.sqrt(Ref_Photon/PhotonsCount)

plt.figure(figsize=(10,8))
plt.subplot(211)
plt.plot(PhotonsCount/Ref_Photon, TransitError, "ko")
plt.xlabel("Normalized Photons Count")
plt.ylabel("Transit Error")
plt.subplot(212)
plt.plot(Wavelength, PhotonsCount/Ref_Photon, "ko")
plt.xlabel("Normalized Photons Count")
plt.ylabel("Transit Error")
plt.tight_layout()
plt.show()


#Now scale the transit error asPhotonsCount/
plt.plot(WavelengthHR, WASP_BlackBody, "k-")
plt.plot(Wavelength, PhotonsMean, "bd")
plt.xlabel("Wavelength [Microns]")
plt.ylabel("Blackbody Flux")
plt.show()



fig, ax = plt.subplots(figsize=(12,8), nrows=2, ncols=1)
plt.subplot(211)
plt.errorbar(Wavelength, TransitDepth, yerr=TransitError, marker="o", capsize=3, linestyle="None")
plt.xlabel("Wavelength [Microns]")
plt.ylabel("Rp/Rs [ppm]")
plt.subplot(212)
plt.errorbar(Wavelength, TransitDepth, yerr=TransitError, marker="o", capsize=3, linestyle="None")
plt.tight_layout()
plt.show()


np.savetxt("WASP_96.txt", np.transpose((WavelengthLower, WavelengthUpper, Wavelength, TransitDepth, TransitError)), header="WavelengthLower, WavelengthUpper, Wavelength, TransitDepth, TransitError")

#Now save the transit 

