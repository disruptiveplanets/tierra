![plot](./figures/LogoLR.png)
Authors: 
1. Prajwal Niraula --- pniraula@mit.edu
2. Aaron Householder --- aaron593@mit.edu 
Insitute: Earth, Atmosphere, and Planetary Science, MIT

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
![plot](./figures/Flowchart.pdf)



