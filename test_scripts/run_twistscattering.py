#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 22:22:23 2020

@author: ramyagurunathan

Twist Scattering Test
"""

import sys
sys.path.append('../scattering_scripts/')

import TwistScattering as TS
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

geom = 'twist'

ax = {'n': 1, 'm' : 2}

d_GS = 350E-09

twist = TS.initialize(input_dict, cmat, density, theta, geom, ax = ax, d_GS = d_GS)

SPlt.convergence_tau_plot(twist, TS.Gamma_rot, 110, T = 300)