#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:50:20 2020

@author: ramyagurunathan

Heterointerface Testing Script
"""

import sys
sys.path.append('../scattering_scripts/')

import HetIntScattering as HS
import ThermalTransport as TT
import ScatteringPlots as SPlt
import numpy as np

input_dict = {
        'avg_vs': (6084 + 5400) / 2, #currently using average velocity between material 1 and material 2
             'atmV': [1.97E-29, 2.27E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

'''
Stiffness matrices for Si and Ge
'''
#Si
cmat1 = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])
density1 = 2330

#Ge
cmat2 = np.array([[126.0, 44.0, 44.0, 0, 0, 0],[44.0, 126.0, 44.0, 0, 0, 0],[44.0, 44.0, 126.0, 0, 0 ,0],\
     [0, 0, 0, 67.7, 0, 0],[0, 0, 0, 0, 67.7, 0],[0 ,0, 0, 0, 0, 67.7]])

density2 = 5323



geom = 'heterointerface'

het = HS.initialize(input_dict, cmat = [cmat1, cmat2], density = [density1, density2], geom = geom)

SPlt.convergence_tau_plot(het, HS.Gamma_rot, 110, T = 300)