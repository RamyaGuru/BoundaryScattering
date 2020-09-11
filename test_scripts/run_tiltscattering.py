#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 21:53:20 2020

@author: ramyagurunathan


Run the Tilt Boundary Scripts
"""

#Specify the tilt boundary angle as a command-line argument
'''
Input dict
'''
import sys
sys.path.append('../scattering_scripts/')

import TiltScattering_update as TS
import ThermalTransport as TT
import ScatteringPlots as SPlt
import numpy as np

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330

theta = 5

geom = 'tilt'

tilt = TS.initialize(input_dict, cmat, density, theta, geom)

SPlt.convergence_tau_plot(tilt, TS.Gamma, 110, T = 300)


