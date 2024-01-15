![plot](./figures/LogoLR.png)
Authors: 
1. Prajwal Niraula --- pniraula@mit.edu
2. Aaron Householder --- aaron593@mit.edu 

Institute: Earth, Atmosphere, and Planetary Science, MIT


## [What is tierra?]
TIERRA (TransmIssion spEctRoscopy of tRansiting plAnets) is a 1D spEctRoscopy code that  adopted and expanded in python by from the original Matlab by Julien de Wit (de Wit et. al 2016, 2018). Currently only 7 molecules were originally included in the code: methane, carbon-monoxide, carbon-dioxide, water, hydrogen, nitrogen and ozone as part of the study Niraula & de Wit et. al 2022). However this has been recently updated to include a wider range of molecules

The default parameters for generating the cross-section are also

The code uses emcee for performing the data.

## [Installation]

First make the clone of the github by:
```
$ git clone https://github.com/disruptiveplanets/tierra
```

To run this project, install it locally using python3. You might need super user privilege to install.

```
$ cd tierra
$ python3 setup.py install
```

## [Updates]

Here are the updates we are working on:
1. Use different PT Profile
2. Add Emission spectroscopy functionality
3. Add better wrappers around fitting
4. Add capabilities of post-processing code.

## [Installation]

We are re-vamping tierra to make it publicly usable for anyone. Currently 


## [Cross-section]

Cross-section should be separately downloaded. All the nine cross-sections related data are around 160GB space in your local space. If you are using the default cross-section this can be directly downloaded from the zenodo at:

The Link will be updated soon



## [Framework]

The framework for perturbation test and is further discussed under the paper Niraula et. al 2022
![plot](./figures/Flowchart.png)


#####################################################################################################################################
For Niraula & de Wit et al. 2022 paper the following grid is used:
Temperature: 100, 110, 120, 130, 140, 160, 180, 200, 230, 260, 290, 330, 370, 410, 460, 510, 580, 650, 730, 810
Log10 Pressure (atm): -5.00, -4.20, -3.80, -3.50, -3.25, -3.05, -2.95, -2.85, -2.75, -2.65, -2.55, -2.45, -2.30, -2.15, -2.00, -1.85,
                      -1.70, -1.55, -1.4, -1.25, -1.1, -0.95, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, 0.10, 0.00, 0.10, 0.20,
                      0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00
#####################################################################################################################################


#####################################################################################################################################
For WASP-39 b retrieval following grid is used:
Temperature: 100, 110, 120, 130, 140, 160, 180, 200, 230, 260, 290, 330, 370, 410, 460, 510, 580, 650, 730, 810, 900, 1000,
             1100, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000
Pressure:-5.00, -4.20, -3.80, -3.50, -3.25, -3.05, -2.95, -2.85, -2.75, -2.65, -2.55, -2.45, -2.30, -2.15, -2.00, -1.85,
         -1.70, -1.55, -1.4, -1.25, -1.1, -0.95, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, 0.10, 0.00, 0.10, 0.20,
         0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00
#####################################################################################################################################


#####################################################################################################################################
For Titan retrievals
Temperature: 70, 84, 98, 112, 126, 140, 154, 168, 182, 196, 210, 224, 238, 252, 266, 280, 294, 308
             322, 336, 350, 364, 378, 392
Pressure:-5.00, -4.20, -3.80, -3.50, -3.25, -3.05, -2.95, -2.85, -2.75, -2.65, -2.55, -2.45, -2.30, -2.15, -2.00, -1.85,
         -1.70, -1.55, -1.4, -1.25, -1.1, -0.95, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, 0.10, 0.00, 0.10, 0.20,
         0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00
#####################################################################################################################################







