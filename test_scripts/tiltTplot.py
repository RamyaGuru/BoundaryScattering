#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:03:19 2020

@author: ramyagurunathan

Tilt Boundary Update:
    Calculate and Plot new Transport Properties
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/scattering_scripts')

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcdefaults()
mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = '18'
mpl.rcParams['text.usetex'] = True
plt.rcParams['legend.title_fontsize'] = '14'

mpl.rcParams['mathtext.fontset'] = 'custom'

mpl.rcParams['mathtext.bf'] = 'Apple Symbols'

'''
Load relaxation times
'''

tilt1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt1spectral_updatetau.npy')
tilt2 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt2spectral_updatetau.npy')
tilt5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt5spectral_updatetau.npy')
tilt8 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt8spectral_updatetau.npy')
tilt10 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt10spectral_updatetau.npy')
tilt15 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/tilt15spectral_updatetau.npy')

'''
Calculate transport coefficients
'''

color = ['xkcd:black', 'xkcd:plum', 'xkcd:blood red', 'xkcd:dark salmon', 'xkcd:silver']
labels = [r'15$^{\circ}$', r'10$^{\circ}$', r'5$^{\circ}$', r'2$^{\circ}$', r'1$^{\circ}$']

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

tilt = AS(**input_dict, geom = 'tilt', theta = 5, ax = 1, d_GS = 350E-9)

tbc1 = []
tbc2 = []
tbc5 = []
tbc8 = []
tbc10 = []
tbc15 = []

for T in [100, 150, 200, 300]:
    transport1 = TT.transport_coeffs_from_tau(tilt, tilt1[0] / tilt.vs, tilt1[1], T)
    transport2 = TT.transport_coeffs_from_tau(tilt, tilt2[0] / tilt.vs, tilt2[1], T)
    transport5 = TT.transport_coeffs_from_tau(tilt, tilt5[0] / tilt.vs, tilt5[1], T)
    transport8 = TT.transport_coeffs_from_tau(tilt, tilt8[0] / tilt.vs, tilt8[1], T)
    transport10 = TT.transport_coeffs_from_tau(tilt, tilt10[0] / tilt.vs, tilt10[1], T)
    transport15 = TT.transport_coeffs_from_tau(tilt, tilt15[0] / tilt.vs, tilt15[1], T)
    tbc1.append((1 / transport1['TBC']) * 1E9)
    tbc2.append((1 / transport2['TBC']) * 1E9)
    tbc5.append((1 / transport5['TBC']) * 1E9)
    tbc8.append((1 / transport8['TBC']) * 1E9)
    tbc10.append((1 / transport10['TBC']) * 1E9)
    tbc15.append((1 / transport15['TBC']) * 1E9)
plt.figure()  
i = 0
for tbc in [tbc15, tbc10, tbc5, tbc2, tbc1]:
    plt.plot([100, 150, 200, 300], tbc, '-o', color = color[i], label = labels[i])
    i = i+1

plt.legend(prop = {'size': 12}, loc = (1.02,0.2))
plt.xlabel('T(K)')
plt.ylabel(r'$R_K$ (m$^2$K/W)')
plt.savefig('TBC_comp_tilt_update.pdf', bbox_inches = 'tight')


'''
Plot the Read Shockley Energy versus the interfacial thermal resistance
'''

tbc_list = [tbc1[-1], tbc2[-1], tbc5[-1], tbc8[-1], tbc10[-1]]

theta = [1,2,5,8,10]

GBenergy = []

for t in theta:
    as_obj = AS(**input_dict, geom = 'tilt', theta = t, ax = 1, d_GS = 350E-9)
    GBenergy.append(AS.gb_energy(as_obj))

plt.figure()
plt.scatter(tbc_list, GBenergy)