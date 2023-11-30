import matplotlib.pyplot as plt
import numpy as np
import glob
import emcee
from tierra import Target, transmission

import os
import sys

import warnings
warnings.filterwarnings("ignore")



def ParseInitiateFile(Location):
    '''
    Give the location of the MCMC initiation files

    Returns:
        The number of parameters in the files given by the number of the lines.
    '''

    FileContent = open(Location).readlines()
    NumParameters = len(FileContent)
    ParametersDict = {}
    for Entry in FileContent:
        #split the values from the entry
        Key, ValueEntries = Entry.split("#")[0].split(":")

        #Check if there is comma:
        if "," in ValueEntries:
            Value, Error = ValueEntries.split(",")
            Value = np.float(Value)
            Error = np.float(Error)
            ParametersDict[Key]= [Value,Error]
        else:
            Value = np.float(ValueEntries)
            ParametersDict[Key] = [Value]

    return NumParameters, ParametersDict



def MoleculesRequired(ParamsDict):
    '''
    Returns the molecules required after parsing the mcmc.ini file.
    '''
    MoleculesExpected = []
    for Key, _ in ParamsDict.items():
        if "MR" in Key:
            Key = Key.replace("MR_","")
            MoleculesExpected.append(Key)
    return MoleculesExpected


def ConvertDict(ParamsDict):
    '''
    Returns the molecules required after parsing the mcmc.ini file.
    '''
    PlanetaryDict = {}
    RequiredParam = ["Mass", "Radius", "T0", "P0", "HeH2Ratio"]
    for key, Value in ParamsDict.items():
        if "MR" in key[:3]:
            PlanetaryDict[key]=Value[0]
        elif "MASS" in key.upper():
            PlanetaryDict['Mass'] = Value[0]
        elif "RADIUS" in key.upper():
            PlanetaryDict['Radius'] = Value[0]
        elif "PT" in key.upper():
            PlanetaryDict['PT'] = Value[0]
        elif "P0" in key.upper():
            PlanetaryDict['P0'] = Value[0]
        elif "T0" in key.upper():
            PlanetaryDict['T0'] = Value[0]
        elif "TINF" in key.upper():
            PlanetaryDict['Tinf'] = Value[0]
        elif "LOGALR" in key.upper():
            PlanetaryDict['ALR'] = 10**Value[0]
        elif "HEH2RATIO" in key.upper():
            PlanetaryDict['HeH2Ratio'] = Value[0]
    #Make sure the required keys in the dictionaries are there
    return PlanetaryDict

#start the walkers based on the file
NumParams, MCMCDict = ParseInitiateFile("mcmc.ini")
Molecules2Look = MoleculesRequired(MCMCDict)
print("The molecules to look are given by:", Molecules2Look)


CurrentLocation = os.getcwd()
if "gridsan" in CurrentLocation:
    print("Running in MIT Supercloud")
    BaseLocation =  "/home/gridsan/ahouseholder/CrossSectionTIERRA"
else:
    print("Running in the Desktop Aekalavya")
    BaseLocation =  "/media/prajwal/LaCie/CrossSectionFromSuperCloud4JWST/CrossSectionTIERRA"



#Look for all the files
PlanetaryDict = ConvertDict(MCMCDict)
print("The planetary dictionary is given by:", PlanetaryDict)


CS_Case = sys.argv[1]
Instr = sys.argv[2]
NSteps = int(sys.argv[3])
PTCase = int(sys.argv[4])

#Use empty for the case of Titan.
StellarDict = {}

Skiprow=1
StellarDict['Mass']=0.93
StellarDict['MassErr']=0.09
StellarDict['Radius']=0.895
StellarDict['RadiusErr']=0.023

PlanetaryDict['Mass']=0.28*317.907
PlanetaryDict['MassErr']=0.03*317.907
PlanetaryDict['Radius']=1.27*11.2089
PlanetaryDict['RadiusErr']=0.04*11.2089

if Instr.upper().replace(" ","") == "NIRCAM":
    print("Loading NIRCAM data.")
    LCFile = "data/WASP39b_NIRCAM.txt" 
elif Instr.upper().replace(" ","") == "NIRISS":
    print("Loading NIRISS data.")
    LCFile = "data/WASP39b_NIRISS.txt"
elif Instr.upper().replace(" ","") == "NIRSPEC_G395H":
    print("Loading NIRSPEC G395H data.")
    LCFile = "data/WASP39b_NIRSPEC_G395H.txt"
elif Instr.upper().replace(" ","") == "PRISM":
    print("Loading NIRSPEC PRISM data.")
    LCFile = "data/WASP39b_PRISM.txt"    
else:
    assert 1==2, "File Not found"

SaveName = "WASP39_MR1_CS"+CS_Case+"_"+Instr+"_PT"+str(PTCase)


if PTCase==1:
    print("Now running PT case of 1")
    from tierra.samplerNumDensity_Isothermal import RunMCMC
elif PTCase==2:
    print("Now running PT case of 2")
    from tierra.samplerNumDensity_Parametric import RunMCMC  



Molecules2Look = sorted(Molecules2Look)
RunMCMC(Molecules2Look, PlanetaryDict, StellarDict, LCFile=LCFile, CSLocation=BaseLocation, AssignedzStep=0.25,  SaveName=SaveName, NSteps=NSteps, CS_Case=CS_Case, NewStart=False)
