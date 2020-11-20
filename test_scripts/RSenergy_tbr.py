#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 14:12:18 2020

@author: ramyagurunathan

Read-Shockley Energy versus the Thermal Boundary Resistance
"""
from ArrayScattering import ArrayScattering as AS
import numpy as np
from math import pi as pi
import matplotlib.pyplot as plt
import ThermalTransport as TT

'''
Materials Parameters
'''
input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = AS(**input_dict, geom = 'twist', theta = 5, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

def core_energy(Ec, b, theta):
    theta = theta * (pi / 180)
    return 2 * np.sin(2 * theta) * (Ec / b)

def strain_energy(Es, b, theta):
    theta = theta * (pi / 180)
    return (Es / b) * np.log(np.sin(2 * theta))


core = []
strain = []
tot = []

Ec = 1
Es = 1
for t in [1, 2, 3, 5, 7, 8, 9, 10, 15]:
    c = core_energy(Ec, twist.b, t)
    s = strain_energy(Es, twist.b, t)
    core.append(c)
    strain.append(s)
    tot.append(c + s)


twist1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist1spectral_updatetau.npy')
twist2 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist2spectral_updatetau.npy')
twist3 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist3spectral_updatetau.npy')
twist5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist5spectral_updatetau.npy')
twist7 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist7spectral_updatetau.npy')
twist8 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist8spectral_updatetau.npy')
twist10 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist10spectral_updatetau.npy')
twist15 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist15spectral_updatetau.npy')

tbc1 = []
tbc2 = []
tbc3 = []
tbc5 = []
tbc7 = []
tbc8 = []
tbc10 = []
tbc15 = []
for T in [100, 150,200,250,300, 500]:
    transport1 = TT.transport_coeffs_from_tau(twist, twist1[0] / twist.vs, twist1[1], T)
    transport2 = TT.transport_coeffs_from_tau(twist, twist2[0] / twist.vs, twist2[1], T)
    transport3 = TT.transport_coeffs_from_tau(twist, twist3[0] / twist.vs, twist3[1], T)
    transport5 = TT.transport_coeffs_from_tau(twist, twist5[0] / twist.vs, twist5[1], T)
    transport7 = TT.transport_coeffs_from_tau(twist, twist7[0] / twist.vs, twist7[1], T)
    transport8 = TT.transport_coeffs_from_tau(twist, twist8[0] / twist.vs, twist8[1], T)
    transport10 = TT.transport_coeffs_from_tau(twist, twist10[0] / twist.vs, twist10[1], T)
    transport15 = TT.transport_coeffs_from_tau(twist, twist15[0] / twist.vs, twist15[1], T)
    tbc1.append((1 / transport1['TBC']) * 1E9)
    tbc2.append((1 / transport2['TBC']) * 1E9)
    tbc3.append((1 / transport3['TBC']) * 1E9)
    tbc5.append((1 / transport5['TBC']) * 1E9)
    tbc7.append((1 / transport7['TBC']) * 1E9)
    tbc8.append((1 / transport8['TBC']) * 1E9)
    tbc10.append((1 / transport10['TBC']) * 1E9)
    tbc15.append((1 / transport15['TBC']) * 1E9)
    
'''
Plot 
'''

fig, ax = plt.subplots()
ax2 = ax.twinx()


