from socket import inet_pton
import numpy as np
import matplotlib.pyplot as plt
import emcee
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from tierra import Target
from tierra.transmission import TransmissionSpectroscopy

from bisect import bisect
import multiprocessing as mp
import os
from scipy.interpolate import interp1d

import glob
import datetime
import h5py
import time


def logLikelihood(theta, Wavelength, Depth, DepthErr):
    '''
    The log likelihood for calculation
    '''
   
    global LeastResidual, ParameterNames, CurrentSaveName, PlanetParamName, \
           nWalkers, StellarParamsDictionary, CurrentSystem
    
    CurrentM, LogP0, T0, V0, HeH2Ratio = theta[:5]
    LogMR = theta[5:]

    P0 = 10**LogP0
    
    PrintFlag = False

    ### Delete the following variables
    if P0<0 or P0>100.:
        print("The case of P0")
        return -np.inf       

    if not(100.<T0<3000.):
        print("The case of T0")
        return -np.inf            
    
    if np.log10(HeH2Ratio)<-3 or np.log10(HeH2Ratio) > 3:
        print("The case for HeH2ratio")
        return -np.inf

    if np.log10(P0)<-10 or np.log10(P0)>2.0:
        print("The case for the base pressure")
        return -np.inf

    if max(LogMR)>-1:
        print("The maximum case for the mixing ratio")
        return -np.inf

    if min(LogMR)<-12.0:
        print("The minimum case of the mixing ratio")
        return -np.inf

    if abs(V0)>10000.0:
        print("The case for the velocity is here...")

    #Calculating the mixing ratio
    MR = 10**LogMR
    MR_H2 = (1.0 - np.sum(MR)) / (1.0 + HeH2Ratio)

    if MR_H2<1e-24 or MR_H2>1.0:
        print("Prior snag on the mixing ratio for hydrogen")
        return -np.inf

    #print("The mixing ratio of hydrogen is given by:", MR_H2)
    #Calculate the PT Profile
    NewPlanetaryDict ={}    
    NewPlanetaryDict['P0'] = P0

    #Arbitrarily assigned
    NewPlanetaryDict['T0'] = T0   
    NewPlanetaryDict['HeH2Ratio'] = HeH2Ratio

    for counter, MoleculeItem in enumerate(PlanetParamName):
        NewPlanetaryDict[MoleculeItem]=MR[counter]
    NewPlanetaryDict['MR_H2'] = MR_H2

    global PlanetParamsDatabase, StellarParamsDatabase
    #remove this later to make things more automatic
    NewPlanetaryDict['Mass'] = CurrentM
    NewPlanetaryDict['Radius'] =  1.25*11.2089 #slightly less than Bonomo+

    CurrentSystem.InitiateSystem(NewPlanetaryDict)
    CurrentSystem.PT_Profile(zStep=0.25, ShowPlot=False)
    T1 = TransmissionSpectroscopy(CurrentSystem, CIA=True)
    T1.CalculateTransmission(CurrentSystem)
    CurrentModel = T1.Spectrum

    if PrintFlag:
        print("The current planetary parameter dictionary is given by: ", NewPlanetaryDict)
        print("The current stellar parameter dictionary is given by: ", StellarParamsDictionary)
        print("The mean molecular mass of the system is: ", CurrentSystem.mu)
        print("The radius of the system is: ", CurrentSystem.Rp)
        print("The mass of the system is: ", CurrentSystem.Mp)
        print("The gravity of the system is:", CurrentSystem.Gp)
        print("Is hydrogen present:", CurrentSystem.H2Present)
        print("The name of the molecules are given by: ", CurrentSystem.MoleculeNames)
        print("And their molecular mass is given by: ", CurrentSystem.MolecularMass)
        print("The number density at the base of the atmosphere is: ", CurrentSystem.nz0)
        print("The mixing ratios for the current system is given by: ", CurrentSystem.MixingRatios)
        input("We will wait here... 100")


    BinnedModel = np.zeros(len(Wavelength))
    global WavelengthLower, WavelengthUpper

    CurrentWavelength = (1.0+V0/299792.458)*CurrentSystem.WavelengthArray*1e4

    counter = 0
    for Wl, Wp in zip(WavelengthLower, WavelengthUpper):
        StartIndex = bisect(CurrentWavelength, Wl)
        StopIndex = bisect(CurrentWavelength, Wp)
        BinnedModel[counter] = np.mean(CurrentModel[StartIndex:StopIndex])
        counter+=1  

    if np.sum(np.isnan(BinnedModel)):
        for x,y in zip(ParameterNames, theta):
            print(x,"  ", y)
        print("There is nan in the binned model.")
        print("The value of theta is given by:", theta)
        return -np.inf

    #Add Gaussian prior in Mass
    #GaussianPriorMass = (CurrentM-PlanetParamsDatabase['Mass'])**2/(PlanetParamsDatabase['MassErr'])**2
    Residual = np.sum(np.power(Depth-BinnedModel,2)/(DepthErr*DepthErr))#+GaussianPriorMass
    ChiSqr = -0.5*Residual

    if Residual<LeastResidual:
       print("The value of Residual is::", Residual)
       print("Saving the best model.")
       LeastResidual = Residual
       with open("MCMCParams/BestParam"+CurrentSaveName+".txt", 'w+') as f:
         f.write("Residual:"+str(Residual)+"\n")
         for key,value in zip(ParameterNames, theta):
              f.write(key+":"+str(value)+"\n")
      
       if 1==1:
          np.savetxt("Figures/BestModel_%s.txt" % CurrentSaveName, np.transpose((Wavelength, BinnedModel)))
          
          import matplotlib.pyplot as plt
          fig, axs = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1]})
          axs[0].errorbar(Wavelength, Depth, yerr=DepthErr, marker=".", markersize=9, fmt='o', color='blue', label='JWST Data', linestyle='', markeredgecolor='k', alpha=0.8)
          axs[0].plot(Wavelength, BinnedModel, "r+-", lw=2, alpha=0.4, label="Model")
          axs[0].set_xlim([min(Wavelength), max(Wavelength)])
          axs[0].set_ylabel("Transit Depth [ppm]", fontsize=18)
          axs[0].legend(loc=1, fontsize = 12)
          axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
          axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
          axs[0].tick_params(which='both', direction='in', top=True, right=True, width=1.5, labelsize=14)
          axs[0].tick_params(which='major', length=6)

          axs[1].set_facecolor('lightgrey')
          axs[1].errorbar(Wavelength, (Depth-BinnedModel)/DepthErr, yerr=np.ones(len(Wavelength)), marker=".", markersize=9, fmt='o', color='blue', label='Data-Model', linestyle='', markeredgecolor='k', alpha=0.8)
          axs[1].axhspan(-1, 1, color='grey', alpha=0.6)
          axs[1].set_xlim([min(Wavelength), max(Wavelength)])
          axs[1].set_xlabel("Wavelength [Microns]", fontsize=18)
          axs[1].set_ylabel("Deviation [$\sigma$]", fontsize=18)
          axs[1].legend(loc=1, fontsize = 12)
          axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
          axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
          axs[1].tick_params(which='both', direction='in', top=True, right=True, width=1.5, labelsize=12)
          axs[1].tick_params(which='major', length=6)

          plt.tight_layout()
          plt.savefig("Figures/CurrentBestModelNew_%s.png" % CurrentSaveName)
          plt.close('all')
          if PrintFlag:
            input("Crash here...")
          print("Best Model Updated. Figure saved.")
    
    return ChiSqr

def RunMCMC(Molecule2Look, PlanetParamsDict, StellarParamsDict, LCFile=None, CSLocation=None, AssignedzStep=0.25,  SaveName="Default", NSteps=50000, NCORES=4, CS_Case="1", NewStart=False, SLICE=True):
    '''
    Run MCMC value.

    Parameters
    ##########

    PlanetParamDict: dictionary
                     Dictionary containing planetary parameter value

    CSLocation: string
                Base location of the cross-section

    AssignedzStep: float
                    Assigned value for the zStep size. Should be smaller than 0.15

    StellarParamDict: dictionary
                      Dictionary containing stellar parameter value

    NumberPTLayers: integer
                    Number of PT layers for the calculation

    NSteps: integer
            Number of steps
    '''


    if "gridsan" in CSLocation:
        print("Using half of the ")
        os.environ["OMP_NUM_THREADS"] = str(mp.cpu_count()//3)
    else:    
        #By default use four cores in CS
        print("Currently overriding this by 4")
        os.environ["OMP_NUM_THREADS"] = "4"#str("NCORES")

    global LeastResidual, CurrentSaveName, CurrentzStep, StartTime, AllMolecules, PlanetParamsDatabase, StellarParamsDatabase

    PlanetParamsDatabase = PlanetParamsDict
    StellarParamsDatabase = StellarParamsDict

    print("Molecules 2 look at line 195", Molecule2Look)
    # This is temporary assignment
    Molecule2Look.append("H2")   
    PlanetParamsDict["MR_H2"]=0.90     
    print("Molecules 2 look at line 199", Molecule2Look)

    AllMolecules = Molecule2Look

    StartTime = str(datetime.datetime.now()).replace(" ", "-").replace(":","-").split(".")[0]
    CurrentSaveName = SaveName
    LeastResidual = np.inf
    CurrentzStep = AssignedzStep
    BaseLocation = CSLocation

    global CurrentSystem, StellarParamsDictionary

    

    StellarParamsDictionary = StellarParamsDict
    PlanetParamsDict['PT']=1
    CurrentSystem = Target.System(Molecule2Look, PlanetParamsDict, StellarParamsDict)

    print("Now loading CS-", CS_Case)
    CurrentSystem.LoadCrossSection(BaseLocation, SubFolder="CS_"+CS_Case, CIA_Flag=True)


    try:
        Wavelength, WavelengthBin, Depth, DepthErr, Order = np.loadtxt(LCFile, unpack=True, skiprows=1, delimiter=",")
    except:
        #This is fire firefly reduction for WASP-39b
        print("Now loading for WASP-39")
        Wavelength, WavelengthBin, Depth, DepthErr = np.loadtxt(LCFile, unpack=True)
        Depth*=1e6
        DepthErr*=1e6
    
    global WavelengthLower, WavelengthUpper
    WavelengthLower = Wavelength-WavelengthBin
    WavelengthUpper = Wavelength+WavelengthBin
   
    #Slice
    if SLICE:
        SmallestWavelength = np.min(WavelengthLower)
        HighestWavelength = np.max(WavelengthUpper)
        StartIndex =  bisect(CurrentSystem.WavelengthArray*1e4, SmallestWavelength)-3
        StopIndex =  bisect(CurrentSystem.WavelengthArray*1e4, HighestWavelength)+300

        CurrentSystem.CrossSectionData = CurrentSystem.CrossSectionData[:,:,StartIndex:StopIndex,:]
        CurrentSystem.WavelengthArray = CurrentSystem.WavelengthArray[StartIndex:StopIndex]

        if hasattr(CurrentSystem, 'CIA_CS'):
            CurrentSystem.CIA_CS = CurrentSystem.CIA_CS[:,:,StartIndex:StopIndex]

    global StartIndexAll, StopIndexAll, nWalkers

   
    #Check for the file name
    FileName = "ProgressData/{}.h5".format(SaveName)
    FileExist = os.path.exists(FileName) #len(FileName)>0 check for the files later 
     
    global backend, sampler
    backend = emcee.backends.HDFBackend(FileName)

    global ParameterNames, PlanetParamName

    #Construct the parta
    ParameterNames = ["Mass", "P0", "T0", "Vplanet", "HeH2Ratio"]
    PlanetParamName = []

    for MoleculeItem in Molecule2Look[:-1]:
        ParameterNames.append("LogMR_"+MoleculeItem)
        PlanetParamName.append("MR_"+MoleculeItem)
    
    print('ParameterNames =', ParameterNames)
    nDim = len(ParameterNames)
    nWalkers = nDim*4

    #Making sure there are at least some saved samples
    if FileExist:
        with h5py.File(FileName, "r") as f:
            LogProb = np.array(f['mcmc']['log_prob'])
            MeanLogProb = np.abs(np.mean(LogProb, axis=1))
            SelectIndex = MeanLogProb>1e-5
            NCases = len(MeanLogProb[SelectIndex])
            print("The number of run steps is:", NCases)
            if NCases<5:
                FileExist=False

    if not(FileExist):
        backend.reset(nWalkers, nDim)

    #Convert the mixing ratio into log of the mixing ratio
    if NewStart or not(FileExist):
        MassInit = np.random.normal(PlanetParamsDict['Mass'],PlanetParamsDict['MassErr'], nWalkers)         #The mass of the exoplanet.
        LogP0Init = np.random.normal(-1.0, 0.5, nWalkers)
        T0Init = np.random.normal(700., 10., nWalkers)
        V0Init = np.random.normal(0., 1000., nWalkers)
        HeH2RatioInit = np.random.uniform(0.1, 0.2, nWalkers)  # Example range for HeH2Ratio

        #Now for molecular number density
        LogMR = np.random.uniform(-3, -7, (nWalkers, nDim-5))     #Mixing ratio for the molecules
        
        StartingGuess = np.column_stack((MassInit, LogP0Init, T0Init, V0Init, HeH2RatioInit, LogMR))
        print("StartingGuess shape:", StartingGuess.shape, "Number of Walkers:", nWalkers, "Number of Dimensions:", nDim)
        print("StartingGuess:", StartingGuess)
    sampler = emcee.EnsembleSampler(nWalkers, nDim, logLikelihood, args=[Wavelength, Depth, DepthErr], backend=backend)
    
    if not(FileExist):
        print("Starting ANEW")
        sampler.run_mcmc(StartingGuess, NSteps, store=True)

    else:
        print("Starting NOT fresh")
        sampler.run_mcmc(None, NSteps, store=True)

    '''with mp.Pool() as pool: 
        sampler = emcee.EnsembleSampler(nWalkers, nDim, logLikelihood, args=[Wavelength, Depth, DepthErr], backend=backend, pool=pool)

        if not(FileExist):
            print("Starting ANEW")
            sampler.run_mcmc(StartingGuess, NSteps, store=True)

        else:
            print("Starting NOT fresh")
            #input("Now starting MCMC")
            sampler.run_mcmc(None, NSteps, store=True)
    '''



