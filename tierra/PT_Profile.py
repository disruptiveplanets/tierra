import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import time
from .Values import Constants

def matlab_multiply(Array_A, Array_B):
    '''perform matlab style matrix multiplication with
    implicit expansion of dimesion'''
    assert np.ndim(Array_A)==1 & np.ndim(Array_B)==1
    NewArray = np.zeros((len(Array_B), len(Array_A)))
    for i in range(len(Array_B)):
        NewArray[i,:] = Array_B[i]*Array_A
    return NewArray


def Analytical_PT_Profile_Corrected(ParamPlanet,  NSpacing=500, Plot=False):
    '''
    This function takes in the parameters of the planet and generates
    the pressure-temperature profile given the temperature and pressure
    specified. In the future, we should be able to take the input from
    GCM.

    ParamPlanet is the planet parameters.

    Plot is the toggle for generating PT profile.
    '''

    Const = Constants()

    #Unpacking the values
    Mp = ParamPlanet['Mass']*Const.M_ear                     #Mass in kg
    Rp = ParamPlanet['Radius'] *Const.R_ear/1000             #Radius in meters
    P0 = ParamPlanet['Pressure']*Const.P_atm                 #Pressure in atmospheric pressure converted to Pascal
    T0 = ParamPlanet['Temperature']                          #Temperature at Rp
    Gam =  ParamPlanet['ALR']                                #Adiabatic Lapse Rate in [K.km^{-1}]
    Tinf =  ParamPlanet['Temperature_Space']                 #Temperature in space in [K]:


    Data = loadmat("DataMatrix/sigma_SuperEarth_with_RS_R_25000.mat")
    Lam = Data['Lam'][0]
    Tem_vect = Data['Tem_vect'][0]
    exp_P = Data['exp_P'][0]

    mixing_ratios = np.array([ParamPlanet['MR_H2O'], ParamPlanet['MR_CO2'],
                              ParamPlanet['MR_O3'], ParamPlanet['MR_CH4'],
                              ParamPlanet['MR_N2'], ParamPlanet['MR_H2']])

    n_z = (1.e-6/Const.k_bo)*mixing_ratios*(P0/T0)           #in cm-3
    n_z[-1] = n_z[-1]*1.176 #To account for the invisible fraction of atoms of He


    #Load the data matrix from
    molparam = loadmat('DataMatrix/molparam.mat')['molparam']
    molparam = np.array(molparam)


    molparam = AddData(molparam)
    list_comp = np.array([1,2,3,6,22,43])-1

    molparam_values = molparam[list_comp,0,0]

    #Need to read the data from molparam
    mu = (sum(molparam_values*n_z)+0.3*n_z[-1])/sum(n_z) # the 0.3 is for the He associated with the H2

    #Is there an analytical solution
    H0 =1.0e6*Const.k_bo*T0/Const.G_gr*(Rp)**2.0/(mu/Const.N_av*Mp)   #taking into account g is km and mu is in grams)

    z_Values = np.linspace(0,25,NSpacing)*H0
    z_Step = z_Values[1] - z_Values[0]
    T_z_analytical = Tinf-(T0-Tinf)*np.exp(-z_Values*Gam)


    #Directly Caluculate the P_z
    g = Mp*Const.G_gr/Rp**2

    #Scale height in kilometers
    Hinf = 1.0e6*Const.k_bo*Tinf/(mu/Const.N_av*g)

    Integral = z_Values/Hinf  \
               - 1./(Hinf*Gam)*np.log(np.abs((Tinf)/(T0-Tinf)+1)) \
               +1./(Hinf*Gam)*np.log(np.abs(Tinf/(T0-Tinf)+np.exp((Gam*z_Values))))

    #Integral = z_Values/Hinf  \
    #           - (T0-Tinf)/(Hinf*Gam)*np.log(((Tinf)/(T0-Tinf)+1)) \
    #           +(T0-Tinf)/(Hinf*Gam)*np.log((Tinf/(T0-Tinf)+np.exp((Gam*z_Values)/(Tinf - T0))))

    P_z_analytical = P0/Const.P_atm*np.exp(-Integral)

    #Select the range of P
    SelectIndex = np.logical_and(P_z_analytical>1e-8,P_z_analytical<1e1)


    #Slice the index
    z_Values = z_Values[SelectIndex]
    T_z_analytical = T_z_analytical[SelectIndex]
    P_z_analytical = P_z_analytical[SelectIndex]

    #coef_P_n_z = P_z_analytical/P0

    #n_z = matlab_multiply(np.array(coef_P_n_z),np.array(n_z))
    #n_z[-1,:] /= 1.176                                                 # Back to H2 only

    if Plot:
        #Generating the figure
        fig, ax = plt.subplots(figsize=(14,6),nrows=1, ncols=2)
        ax[0].plot(P_z_analytical,z_Values, "r-", linewidth=1)
        ax[0].set_xlabel('Pressure (atm)', color='red', fontsize=20)
        ax[0].set_ylabel('Atmosphere (km)', color='blue', fontsize=20)
        ax[0].grid(True)
        ax[0].tick_params(axis='x', labelcolor='red')
        ax[0].set_xscale('log')

        ax_0 = ax[0].twiny()
        ax_0.plot(T_z_analytical,z_Values, "g-", linewidth=1)
        ax_0.set_xlabel('Temperature (K)', color='green', fontsize=20)
        ax_0.tick_params(axis='x', labelcolor='green')


        #Name of the molecules
        LineStyles = [':', '-.', '--', '-']
        for i in range(np.shape(n_z)[0]):
            ax[1].plot(n_z[i,:],z_Values,linewidth=1, linestyle=LineStyles[i%len(LineStyles)], label=MolName[i])
        ax[1].grid(True)
        ax[1].set_xscale('log')
        ax[1].set_ylabel('Atmosphere (km)', color='blue', fontsize=20)
        ax[1].legend(loc=1)
        plt.tight_layout()
        plt.show()

    return [z_Values, n_z, T_z_analytical, P_z_analytical]
