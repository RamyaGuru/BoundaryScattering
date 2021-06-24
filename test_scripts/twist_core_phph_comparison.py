#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:26:47 2021

@author: ramyagurunathan

Twist boundary: comparison to model with core and phonon-phonon interactions
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

twist = AS(**input_dict, geom = 'twist', theta = 11, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

twist3 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/lit_angles/twist3.4spectral_update2tau.npy') 
twist3_phph = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist3.4spectral_phphtau.npy')
twist3_phph_core = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist3.4spectral_phph_coretau.npy')



plt.plot(twist3[0] / twist.omegaD, twist3[1], label = 'strain only')
plt.plot(twist3_phph[0] / (twist.omegaD * (pi / 2)), twist3_phph[1], label = 'strain + phph')
plt.plot(twist3_phph_core[0] / (twist.omegaD * (pi / 2)), twist3_phph_core[1], label = 'strain + phph + core')
plt.legend()


'''
Loglog plots
'''

plt.figure()
#plt.loglog(twist3[0] / twist.omegaD, twist3[1], label = 'strain only')
plt.loglog(twist3_phph[0] / (twist.omegaD * (pi / 2)), twist3_phph[1], label = 'strain + phph')
plt.loglog(twist3_phph[0][4:7] / (twist.omegaD * (pi / 2)), twist3_phph[1][4]*(twist3_phph[0][4:7]/twist3_phph[0][4])**(-4) * 0.75, color = 'xkcd:black')
plt.loglog(twist3_phph[0][50:] / (twist.omegaD * (pi / 2)), twist3_phph[1][50]*(twist3_phph[0][50:]/twist3_phph[0][50])**(-1.1) * 0.75, color = 'xkcd:black')
plt.loglog(twist3_phph_core[0] / (twist.omegaD * (pi / 2)), twist3_phph_core[1], label = 'strain + phph + core')
plt.loglog(twist3_phph_core[0][50:] / (twist.omegaD  * (pi / 2)), twist3_phph_core[1][50]*(twist3_phph_core[0][50:]/twist3_phph_core[0][50])**(-3) * 0.75, color = 'xkcd:black')
plt.loglog(twist3_phph_core[0][3:7] / (twist.omegaD  * (pi / 2)), twist3_phph_core[1][3]*(twist3_phph_core[0][3:7]/twist3_phph_core[0][3])**(-9) * 0.75, color = 'xkcd:black')
plt.legend()


twist7_phph_core = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist7.0spectral_phph_coretau.npy')
plt.figure()
plt.plot(twist7_phph_core[0] / (twist.omegaD * (pi / 2)), twist7_phph_core[1])
plt.figure()
plt.loglog(twist7_phph_core[0] / (twist.omegaD * (pi / 2)), twist7_phph_core[1])

plt.loglog(twist7_phph_core[0][50:] / (twist.omegaD  * (pi / 2)), twist7_phph_core[1][50]*(twist7_phph_core[0][50:]/twist7_phph_core[0][50])**(-1.1)* 0.9, color = 'xkcd:black')
plt.loglog(twist7_phph_core[0][9:14] / (twist.omegaD  * (pi / 2)), twist7_phph_core[1][9]*(twist7_phph_core[0][9:14]/twist7_phph_core[0][9])**(-3.1), color = 'xkcd:black')

plt.figure()
plt.plot(twist7_phph_core[0] / (twist.omegaD * (pi / 2)), twist7_phph_core[1])
plt.xlim([0, 0.1])
plt.ylim([3.5, 4.0])

plt.figure()
plt.loglog(twist7_phph_core[0][0:8] / (twist.omegaD  * (pi / 2)), twist7_phph_core[1][0:8])
plt.loglog(twist7_phph_core[0][0:8] / (twist.omegaD  * (pi / 2)), twist7_phph_core[1][0]*(twist7_phph_core[0][0:8]/twist7_phph_core[0][0])**(-1), color = 'xkcd:black')



twist7_core_only = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist7.0core_onlytau.npy')

plt.figure()
plt.plot(twist7_core_only[0] / (twist.omegaD), twist7_core_only[1])  #not doing much...

plt.figure()
plt.loglog(twist7_core_only[0] / (twist.omegaD * (pi / 2)), twist7_core_only[1])
#plt.loglog(twist7_core_only[0][80:] / (twist.omegaD  * (pi / 2)), twist7_core_only[1][80]*(twist7_core_only[0][80:]/twist7_core_only[0][80])**(-0.01), color = 'xkcd:black')
