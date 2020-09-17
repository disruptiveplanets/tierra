import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect
import sys

class TransmissionSpectroscopy:

    def __init__(self, Target):
        '''
        Initiate the transmission
        '''



        '''
        self.dz = np.diff(Target.zValues)
        Z_ii, Z_jj = np.meshgrid(Target.zValues[1:], Target.zValues[:-1])
        self.Distance = np.sqrt((Target.Rp+Z_ii)*(Target.Rp+Z_ii)
                   -(Target.Rp+Z_jj)*(Target.Rp+Z_jj))

        self.Distance[np.isnan(self.Distance)]=0.0
        '''


        sz = Target.NumLayers
        z_ = Target.zValuesCm
        self.dz_= np.concatenate(([z_[0]], z_[1:sz]-z_[:sz-1]))
        x__= np.zeros((sz-1,sz-1))

        for i in range(sz-1):
            for j in range(i,sz-1):
                x__[i,j]=np.sqrt((Target.Rp+Target.zValuesCm[j+1])**2.0-(Target.Rp+Target.zValuesCm[i])**2.0)




        XNew = np.zeros(np.shape(x__))
        XNew[:,1:] = x__[:,:sz-2]

        self.ds_=x__- XNew
        self.ds_/=1.e5
        #ds_[ds_==0] = np.nan

        XNew_ = np.zeros(np.shape(self.ds_))
        XNew_[:,1:] = self.ds_[:,:sz-2]
        self.ds_ = 0.5*(self.ds_ + XNew_)


        #Initiating the alpha function
        self.alpha = np.zeros((len(Target.WavelengthArray),Target.NumLayers),dtype=np.float32)





    def CalculateTransmission(self, Target, ShowPlot=False):
        '''
        This method calculates the spectrum given the

        Parameters:
        -----------
        Target: Tierra Target object


        ShowPlot: Boolean

        Returns
        --------
        Array

        Spectrum of the planet is returned.

        '''

        #Now solve for the atmosphere of the planet

        self.Spectrum = np.zeros(len(Target.WavelengthArray), dtype=np.float32)



        for self.CurrentLayer in range(Target.NumLayers):
            CurrentT = Target.TzAnalytical[self.CurrentLayer]
            CurrentP = np.log10(Target.PzAnalytical[self.CurrentLayer])

            TIndex = bisect(Target.TemperatureArray, CurrentT)
            PIndex = bisect(Target.PressureArray, CurrentP)

            Temp1, Temp2 = [Target.TemperatureArray[TIndex-1], Target.TemperatureArray[TIndex]]
            P1, P2 = [Target.PressureArray[PIndex-1], Target.PressureArray[PIndex]]

            #See if they are exactly same
            m = (CurrentP-P1)/(P2-P1)

            #Assign 0
            if not(Target.SmallFile):
                self.alpha[:,self.CurrentLayer] = 0.0
                for Counter, self.CurrentMolecule in enumerate(Target.MoleculeName):
                    print("The name of the molecule is ::", self.CurrentMolecule)

                    if "CH4" in self.CurrentMolecule:
                        CurrentCS = Target.CH4Data
                    elif "CO2" in self.CurrentMolecule:
                        CurrentCS = Target.CO2Data
                    elif "CO" in self.CurrentMolecule:
                        CurrentCS = Target.COData
                    elif "H2O" in self.CurrentMolecule:
                        CurrentCS = Target.H2OData
                    elif "H2" in self.CurrentMolecule:
                        CurrentCS = Target.H2Data
                    elif "O3" in self.CurrentMolecule:
                        CurrentCS = Target.O3Data
                    elif "N2" in self.CurrentMolecule:
                        CurrentCS = Target.N2Data
                    else:
                        print("No cross-section found...")
                        continue


                    Sigma11 = CurrentCS[TIndex-1, PIndex-1, :]
                    Sigma12 = CurrentCS[TIndex-1, PIndex, :]
                    Sigma21 = CurrentCS[TIndex, PIndex-1, :]
                    Sigma22 = CurrentCS[TIndex, PIndex, :]


                    #Performing hill interpolation
                    UndSigma1 = Sigma11+ m*(Sigma12-Sigma11)
                    UndSigma2 = Sigma21+ m*(Sigma22-Sigma21)

                    RatioSigma = UndSigma1/UndSigma2
                    bi = 1./(1./Temp2-1./Temp1)*np.log(RatioSigma)
                    ai = UndSigma1*np.exp(bi/Temp1)

                    self.CurrentInterpSigma = ai*np.exp(-bi/CurrentT)

                    ZeroIndex = np.logical_or(np.isnan(self.CurrentInterpSigma), \
                                              ~np.isfinite(self.CurrentInterpSigma))

                    print("Interpolated Index::")
                    print("Invalid values:", np.sum(self.CurrentInterpSigma>1e-5))

                    print("The current Interpolated sig")
                    #Replace nan with zeros
                    self.CurrentInterpSigma[ZeroIndex] = 0.0

                    #Multiply the cross section with the molecular number density
                    self.alpha[:,self.CurrentLayer] += self.CurrentInterpSigma*Target.nz[Counter,self.CurrentLayer]

                    #self.alpha[:,self.CurrentLayer] += Sigma11*Target.nz[Counter,self.CurrentLayer]


                    #The following is just for some diagnostic plot...
                    if CurrentP>-5.0:

                        Sampling = 5
                        Wavelength = Target.WavelengthArray*1e7
                        SelectIndex = Wavelength<5000


                        Y1 = self.CurrentInterpSigma
                        Y2 = self.alpha[:,self.CurrentLayer]

                        plt.figure(figsize=(14,12))
                        plt.subplot(311)
                        plt.plot(Wavelength[SelectIndex], Y1[SelectIndex] , "k-", lw=0.5, label="Interpolated")
                        plt.title(self.CurrentMolecule)
                        plt.subplot(312)
                        plt.plot(Wavelength[SelectIndex], Sigma11[SelectIndex], "g-", lw=0.5, label="Interpolated")
                        plt.subplot(313)
                        plt.plot(Wavelength[SelectIndex], Y2[SelectIndex], "r-", lw=0.5, label="Alpha")
                        plt.legend()
                        plt.tight_layout()
                        SaveName = "DiagnosticPlots/Layer"+str(self.CurrentLayer).zfill(3)+"_Mol" \
                                    +str(Counter).zfill(3)+"_"+str(int(CurrentT))+"K"+\
                                    str(round(CurrentP,5))+"_interp.png"
                        plt.savefig(SaveName)
                        plt.close('all')

            elif Target.SmallFile:

                co_t = (CurrentT-Target.TemperatureArray[TIndex-1])/(Target.TemperatureArray[TIndex]-Target.TemperatureArray[TIndex-1])

                #co_p = (10.**CurrentP-10.**Target.PressureArray[PIndex-1])/(10.**Target.PressureArray[PIndex]-10.**Target.PressureArray[PIndex-1])

                co_p = (CurrentP-Target.PressureArray[PIndex-1])/(Target.PressureArray[PIndex]-Target.PressureArray[PIndex-1])


                FirstTerm = np.matmul(Target.CrossSectionData[TIndex-1, PIndex-1,:,:], Target.nz[:, self.CurrentLayer])
                SecondTerm = np.matmul(Target.CrossSectionData[TIndex-1, PIndex,:,:], Target.nz[:, self.CurrentLayer])
                ThirdTerm = np.matmul(Target.CrossSectionData[TIndex, PIndex-1,:,:], Target.nz[:, self.CurrentLayer])
                FourthTerm = np.matmul(Target.CrossSectionData[TIndex, PIndex,:,:], Target.nz[:, self.CurrentLayer])

                self.alpha[:,self.CurrentLayer] = ((1-co_t)*(1-co_p))*FirstTerm + \
                                                  ((1-co_t)*co_p)*SecondTerm +  \
                                                  (co_t*(1-co_p))*ThirdTerm + \
                                                  (co_t*co_p)*FourthTerm





        R_s = 1.0e5
        sz = Target.NumLayers
        self.Spectrum = ((Target.Rp+ Target.zValuesCm[0])**2+ \
                        2.0*np.matmul(1.0-(np.exp(-(2.0*(np.matmul(self.alpha[:,0:sz-1],np.transpose(self.ds_[:,:sz-1])))))), \
                        (Target.Rp+Target.zValuesCm[:sz-1])*np.transpose(self.dz_[:sz-1])))/R_s**2
        self.Spectrum = self.Spectrum.flatten()
        Target.Rp/=1.e5
        self.SpectrumHeight = 0.5*(self.Spectrum/Target.Rp-Target.Rp)
        return self.SpectrumHeight
