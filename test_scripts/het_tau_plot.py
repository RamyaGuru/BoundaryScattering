#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:59:09 2020

@author: ramyagurunathan

Plot heterointerface results
"""

import numpy as np
import matplotlib.pyplot as plt

#Heterointerface Results

het_tau = np.load('../datafiles/fall2020/heterointerfaceNonespectral_updatetau.npy')

plt.plot(het_tau[0,:], het_tau[1,:])


#Tilt 5 Boundary Results
plt.figure()
tilt5_tau = np.load('../datafiles/fall2020/tilt5spectral_updatetau.npy')
plt.plot(tilt5_tau[0,:], tilt5_tau[1,:])

#Tilt 1 Boundary Results
plt.figure()
tilt1_tau = np.load('../datafiles/fall2020/tilt1spectral_updatetau.npy')
plt.plot(tilt1_tau[0,:], tilt1_tau[1,:])

#Tilt 10 Boundary Results
plt.figure()
tilt10_tau = np.load('../datafiles/fall2020/tilt10spectral_updatetau.npy')
plt.plot(tilt10_tau[0,:], tilt10_tau[1,:])