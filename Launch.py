import matplotlib.pyplot as plt

from tierra import Target
from tierra.transmission import TransmissionSpectroscopy
from tierra.sampler import EmceeSampler

Target = Target.System()

Target.LoadCrossSection("/media/prajwal/a66433b1-e5b2-467e-8ebf-5857f498dfce/Mar17_2020/", SubFolder="CS_1")
#Target.LoadCrossSection("/media/prajwal/a66433b1-e5b2-467e-8ebf-5857f498dfce/LowerResolutionData/R1000", SubFolder="CS_1")
Target.PT_Profile(NumLayers=20, ShowPlot=False)

T1 = TransmissionSpectroscopy(Target)

print("Now Calculating the transmission spectra")
T1.CalculateTransmission(Target)


fig, ax = plt.subplots(figsize=(14,10), nrows=2, ncols=1, sharex=True)
ax[0].plot(Target.WavelengthArray*1e7, T1.Spectrum, "k-")
ax[0].set_xlim([300,5000])
ax[0].set_xscale('log')
ax[0].set_xlabel("Wavelength (nm)")
ax[0].set_ylabel("Spectrum")

ax[1].plot(Target.WavelengthArray*1e7, T1.SpectrumHeight, "r-")
ax[1].set_xlim([300,5000])
ax[1].set_xlabel("Wavelength (nm)")
ax[1].set_ylabel("Height")
ax[1].set_xscale('log')

plt.tight_layout()
plt.savefig("DeleteMe.png")
plt.close()
