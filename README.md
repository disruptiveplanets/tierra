![plot](./logo/LogoLR.png)
Author: Prajwal Niraula
Email: pniraula@mit.edu
Insitute: MIT

## [What is tierra?]
TIERRA (TransmIssion spEctRoscopy of tRansiting plAnets) is a 1D spEctRoscopy code that  adopted and expanded in python by from the original Matlab by Julien de Wit (de Wit et. al 2016, 2018). Currently only 7 molecules are included in the code: methane, carbon-monoxide, carbon-dioxide, water, hydrogen, nitrogen and ozone as part of the study Niraula & de Wit et. al 2022).

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


## [Cross-section]

Cross-section should be separately downloaded. All the nine cross-sections related data are around 160GB space in your local space. If you are using the default cross-section this can be directly downloaded from the zenodo at:
