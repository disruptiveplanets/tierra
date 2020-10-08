import numpy as np
import matplotlib.pyplot as plt
import os
import glob


def matlab_multiply(Array_A, Array_B):
    '''perform matlab style matrix multiplication with
    implicit expansion of dimesion'''
    assert np.ndim(Array_A)==1 & np.ndim(Array_B)==1
    NewArray = np.zeros((len(Array_B), len(Array_A)))
    for i in range(len(Array_B)):
        NewArray[i,:] = Array_B[i]*Array_A
    return NewArray


class System:
    def __init__(self, PlanetParamsDict=None, StellarParamsDict=None, LoadFromFile=True):
        '''
        LoadFromFile: bool
                      True if to be loaded from the data

        '''

        self.InitiateConstants()
        self.LoadFromFile = LoadFromFile

        self.MolDict = {"H2O":18.010565, "CO2":43.989830, "O3":47.984745,  "N2O":44.001062, \
        "HCl":35.976678, "CO":27.994915,  "CH4":16.031300, "NH3":17.026549, \
        "O2":31.98983,   "H2":2.015650,   "He":4.002602,   "N2":28.006148}

        self.MoleculeName = np.array(['H2O', 'CO2', 'CO', 'O3', 'CH4', 'N2', 'H2'])

        if not(StellarParamsDict):
            print("Assigning stellar parameters from the files.")
            self.ParseStarFile()
        else:
            print("Assigning stellar parameters from the dictionary.")
            self.StellarParams = {}
            for key, value in StellarParamsDict.items():
                self.StellarParams[key] = value

        if self.LoadFromFile:
            self.InitiateSystem()
        elif PlanetParamsDict:
            self.PlanetParams = {}
            for key, value in PlanetParamsDict.items():
                self.PlanetParams[key] = value


    def InitiateSystem(self):
        '''
        Initiate the calculation for mean molecular mass
        and assign the value for pressure and temperature
        '''
        if self.LoadFromFile:
            self.ParsePlanetFile()


        #Unpacking the values
        self.Mp = self.PlanetParams['Mass']*self.M_ear                     #Mass in kg
        self.Rp = self.PlanetParams['Radius'] *self.R_ear                  #Radius in meters
        self.Gp = self.G_gr*self.Mp/(self.Rp*self.Rp)
        self.Ms = self.StellarParams['Mass']*self.M_sun
        self.Rs = self.StellarParams['Radius']*self.R_sun

        self.RpKm = self.Rp/1.e5
        self.RsKm = self.Rs/1.e5
        self.P0 = self.PlanetParams['P0']*self.P_atm              #Pressure in atmospheric pressure converted to Pascal
        self.T0 = self.PlanetParams['T0']                         #Temperature at Rp
        self.Gam =  self.PlanetParams['ALR']                      #Adiabatic Lapse Rate in [K.km^{-1}]
        self.Tinf =  self.PlanetParams['TInf']                    #Temperature in space in [K]:


        self.MixingRatios = np.array([self.PlanetParams['MR_H2O'], self.PlanetParams['MR_CO2'], \
                                      self.PlanetParams['MR_CO'], self.PlanetParams['MR_O3'], \
                                      self.PlanetParams['MR_CH4'], self.PlanetParams['MR_N2'], \
                                      self.PlanetParams['MR_H2']])

        self.MolParamValues = np.array([self.MolDict['H2O'], self.MolDict['CO2'],self.MolDict['CO'],\
                                        self.MolDict['O3'], self.MolDict['CH4'], self.MolDict['N2'],\
                                        self.MolDict['H2']])


        #Number density in calculating the mean molecular mass, Take into consideration the invisible Helium fraction
        MuMixingRatio = np.concatenate((self.MixingRatios, [self.MixingRatios[-1]*0.33]))
        MuMolParamValues = np.concatenate((self.MolParamValues, [self.MolDict['He']]))
        MuNumDensity = (1.e-6/self.k_bo)*MuMixingRatio*(self.P0/self.T0)           #in cm-3


        #Mean molecular mass
        self.mu = sum(MuMolParamValues*MuNumDensity)/sum(MuNumDensity)

        #Has the data been loaded.
        self.CSDataLoaded = False


    def ParsePlanetFile(self):
        '''
        This function parses the planetary file
        '''
        self.PlanetParams = {}

        if os.path.exists("PlanetParam.ini"):
            FileContent = open("PlanetParam.ini", "r").readlines()
            for Line in FileContent:
                Item = Line.split("#")[0].replace(" ","")
                key, Value = Item.split(":")
                self.PlanetParams[key]=float(Value)
        else:
            print("PlanetParam.ini does not exist in the local dictionary")


    def ParseStarFile(self):
        '''
        This function parses the star file i.e StelarParam.ini
        '''

        self.StellarParams = {}

        if os.path.exists("StellarParam.ini"):
            FileContent = open("StellarParam.ini", "r").readlines()
            for Line in FileContent:
                Item = Line.split("#")[0].replace(" ","")
                key, Value = Item.split(":")
                self.StellarParams[key]=float(Value)
        else:
            print("StellarParam.ini does nopt exist in the local dictionary")


    def InitiateConstants(self):
        #Astronomical Constants
        self.R_sun=6.957E10                     #Sun radius in cm
        self.T_sun=5770                         #Effective Temperature of the Sun
        self.M_sun=1.98911e33                   #Mass of the sun in kilograms
        self.P_terre=86400*365.25696            #Days in seconds
        self.R_ear=6.371e8                      #earth's radius in centimeters
        self.r_t=1.496e13                       #1 AU in centimeters
        self.parsec = 3.085677e18               #parsec in centimeters
        self.M_jup = 1.8986e30                  #Mass of Jupiter in kilograms
        self.M_ear = 5.9736e27                  #Mass of the earth in kilograms
        self.R_jup = 6.995e9                    #Radius of Jupiter in meters

        #Physical constants
        self.G_gr = 6.67384E-8                  #Gravitational constant in CGS
        self.c_li = 2.998e10                    #Speed of Light in CGS
        self.h_pl = 6.626069E-27                #Planck's Constant in CGS
        self.k_bo = 1.38065E-16                 #Boltzmann Constant in CGS
        self.P_atm= 1.013e6                     #1 atmospheric pressure of earth in CGS
        self.N_av = 6.02214E23                  #Avogadro's Number
        self.sigma_bo = 5.670E-5                #stefan Boltzmann's Constant


    def PT_Profile(self, zStep=0.25, ShowPlot=False):
        '''
        This method calculates the Pressure Temperature profile
        for the planet

        Parameters:
        -----------------

        zStep: Float
                   Stepsize in atmospheric scale.

        PlotFlag: Boolean
                  Default value is False. Plot the data if True.

        '''

        assert np.sign(self.Tinf-self.T0) == -np.sign(self.Gam)

        #atmospheric scaled height in km
        self.H0 =self.k_bo*self.T0/(self.mu/self.N_av*self.Gp)/1e5
        self.zValues = np.arange(0,15,zStep)*self.H0
        self.zValuesCm=self.zValues*1e5
        self.zStep = self.zValues[1] - self.zValues[0]
        self.TzAnalytical = self.Tinf-(self.Tinf-self.T0)*np.exp(self.zValues*self.Gam/(self.Tinf-self.T0))


        #Scale height in kilometers
        self.Hinf = self.k_bo*self.Tinf/(self.mu/self.N_av*self.Gp)/1e5


        self.Integral = self.zValues/self.Hinf  \
                        + (self.T0-self.Tinf)/(self.Hinf*self.Gam)*np.log(np.abs((self.T0)/(self.T0-self.Tinf))) \
                        - (self.T0-self.Tinf)/(self.Hinf*self.Gam)*np.log( np.abs(self.Tinf/(self.T0-self.Tinf) \
                        + np.exp((-self.Gam*self.zValues)/(self.T0 - self.Tinf))))


        self.PzAnalytical = self.P0/self.P_atm*np.exp(-self.Integral)

        SelectIndex = np.logical_and(self.PzAnalytical>1e-8,self.PzAnalytical<1e2)
        #Slice the index
        self.zValues = self.zValues[SelectIndex]
        self.TzAnalytical = self.TzAnalytical[SelectIndex]
        self.PzAnalytical = self.PzAnalytical[SelectIndex]


        #Number density in per cm^3
        self.nz0 = self.N_av*self.P0/self.P_atm*273.15/self.T0*self.MixingRatios
        self.nz0 = self.N_av/22400.0*self.P0/self.P_atm*273.15/self.T0*self.MixingRatios
        self.NumLayers = len(self.zValues)

        self.nz = np.zeros((len(self.nz0), len(self.PzAnalytical)))

        for i in range(len(self.nz0)):
            self.nz[i,:] = self.nz0[i]*self.PzAnalytical

        self.nzSquared = self.nz*self.nz


        if ShowPlot:
            #Generating the figure
            fig, ax = plt.subplots(figsize=(14,6),nrows=1, ncols=2)
            ax[0].plot(self.PzAnalytical, self.zValues, "r-", linewidth=1)
            ax[0].set_xlabel('Pressure (atm)', color='red', fontsize=20)
            ax[0].set_ylabel('Atmosphere (km)', color='blue', fontsize=20)
            ax[0].grid(True)
            ax[0].tick_params(axis='x', labelcolor='red')
            ax[0].set_xscale('log')

            ax_0 = ax[0].twiny()
            ax_0.plot(self.TzAnalytical, self.zValues, "g-", linewidth=1)
            ax_0.set_xlabel('Temperature (K)', color='green', fontsize=20)
            ax_0.tick_params(axis='x', labelcolor='green')


            #Name of the molecules
            LineStyles = [':', '-.', '--', '-']
            for i in range(len(self.nz0)):
                ax[1].plot(self.nz[i,:], self.zValues,linewidth=1, \
                linestyle=LineStyles[i%len(LineStyles)], label=self.MoleculeName[i])
            ax[1].grid(True)
            ax[1].set_xscale('log')
            ax[1].set_ylabel('Atmosphere (km)', color='blue', fontsize=20)
            ax[1].legend(loc=1)
            plt.tight_layout()
            plt.show()



    def LoadCrossSection(self, Location="", SubFolder=""):
        '''
        This method is supposed to load the cross-section

        The expected location is
        '''

        CombinedLocation = os.path.join(Location, SubFolder)
        print("Loading the cross-section from: ", CombinedLocation)

        AllFileList = np.array(glob.glob(os.path.join(Location, SubFolder)+"/*.npy"))
        MoleculeFileList = np.array([FileItem.split("/")[-1][:-4] for FileItem in AllFileList])
        Size = os.path.getsize(os.path.join(CombinedLocation,"H2O.npy"))/1e6

        assert os.path.exists(Location+"/Wavelength.npy"), "Wavelength.npy is needed "
        assert len(MoleculeFileList)>=len(self.MoleculeName), "The number number of molecules are not here."

        self.WavelengthArray = np.load(Location+"/Wavelength.npy")
        self.TemperatureArray = np.loadtxt(Location+"/Temperature.txt")
        self.PressureArray = np.loadtxt(Location+"/Pressure.txt")

        NumWavelength = len(self.WavelengthArray)
        NumTemp = len(self.TemperatureArray)
        NumPressure = len(self.PressureArray)

        print("The size of single molecule ::", Size )
        #If greater than 2 GB
        if Size>4000:
            self.SmallFile = False
        else:
            self.SmallFile = True

        if not(self.SmallFile):
            self.CH4Data = np.load(os.path.join(CombinedLocation, "CH4.npy"), mmap_mode='r')
            self.COData = np.load(os.path.join(CombinedLocation, "CO.npy"), mmap_mode='r')
            self.CO2Data = np.load(os.path.join(CombinedLocation, "CO2.npy"), mmap_mode='r')
            self.H2Data = np.load(os.path.join(CombinedLocation, "H2.npy"), mmap_mode='r')
            self.H2OData = np.load(os.path.join(CombinedLocation, "H2O.npy"), mmap_mode='r')
            self.O3Data = np.load(os.path.join(CombinedLocation, "O3.npy"), mmap_mode='r')
            self.N2Data = np.load(os.path.join(CombinedLocation, "N2.npy"), mmap_mode='r')
            self.CSDataLoaded = True
        elif self.SmallFile:
            self.CrossSectionData = np.zeros((NumTemp, NumPressure, NumWavelength, len(self.MoleculeName)))
            for MolCounter, Molecule in enumerate(self.MoleculeName):
                print(MolCounter, ": ", Molecule)
                CurrentData = np.load(os.path.join(CombinedLocation, Molecule+".npy"), mmap_mode='r')
                self.CrossSectionData[:,:,:,MolCounter] = CurrentData
        pass
