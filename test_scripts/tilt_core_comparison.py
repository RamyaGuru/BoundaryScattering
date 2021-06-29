#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 22:34:44 2021

@author: ramyagurunathan

tilt boundary comparison
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../scattering_scripts/')
from ArrayScattering import ArrayScattering as AS
from math import pi as pi

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'shearmod' : 60E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }


tilt = AS(**input_dict, geom = 'tilt', theta = 11, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

tilt7_core_only = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/tilt7.0core_onlytau.npy')

tilt7_core_strain = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/tilt7core_straintau.npy')

plt.figure()
plt.plot(tilt7_core_only[0] / (tilt.omegaD), tilt7_core_only[1]) 

plt.figure()
plt.loglog(tilt7_core_only[0] / (tilt.omegaD * (pi / 2)), tilt7_core_only[1])
plt.loglog(tilt7_core_only[0][80:] / (tilt.omegaD  * (pi / 2)), tilt7_core_only[1][80]*(tilt7_core_only[0][80:]/tilt7_core_only[0][80])**(-2), color = 'xkcd:black')

plt.figure()
plt.plot(tilt7_core_strain[0] / (tilt.omegaD), tilt7_core_strain[1]) 


'''
tilt boundary: theta = 3
'''
tilt3_core_only = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/tilt3core_onlytau.npy')

tilt3_core_strain = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/tilt3core_straintau.npy')

plt.figure()
plt.plot(tilt3_core_only[0] / (tilt.omegaD), tilt3_core_only[1]) 

plt.figure()
plt.plot(tilt3_core_strain[0] / (tilt.omegaD), tilt3_core_strain[1])




