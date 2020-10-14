#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:57:21 2020

@author: ramyagurunathan

Tilt scattering versus temperature plots
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/scattering_scripts')

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
import matplotlib.ticker as ticker

mpl.rcdefaults()
mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = '16'
mpl.rcParams['text.usetex'] = True
plt.rcParams['legend.title_fontsize'] = '14'

mpl.rcParams['mathtext.fontset'] = 'custom'

mpl.rcParams['mathtext.bf'] = 'Apple Symbols'


'''
Load tilt obundary data
'''

tilt1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/tilt1spectral_updatetau.npy')
tilt5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/tilt5spectral_updatetau.npy')

'''
Load twist boundary data
'''

twist1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist1spectral_updatetau.npy')
twist5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist5spectral_updatetau.npy')

'''
Calculate transport coefficients
'''

labels = [ r'5$^{\circ}$', r'1$^{\circ}$',]

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }


tilt = AS(**input_dict, geom = 'tilt', theta = 5, ax = 1, d_GS = 350E-9, bvK = True)
twist = AS(**input_dict, geom = 'twist', theta = 5, ax = {'n' : 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

tbc1 = []
tbc5 = []

for T in [100, 150,200,250,300]:
    transport1 = TT.transport_coeffs_from_tau(tilt, tilt1[0] / tilt.vs, tilt1[1], T)
    transport5 = TT.transport_coeffs_from_tau(tilt, tilt5[0] / tilt.vs, tilt5[1], T)
    tbc1.append((1 / transport1['TBC']) * 1E9)
    tbc5.append((1 / transport5['TBC']) * 1E9)

plt.figure()

i = 0
for tbc in [tbc5, tbc1]:
    plt.plot([100,150,200,250,300], tbc, '-o', label = labels[i])
    i = i+1
    
plt.legend()

'''
Spectral Tau Plots
'''

plt.rcParams["figure.figsize"] = [5, 3]

mpl.rcParams['font.size'] = '14'
plt.figure()
ax = plt.gca()


tilt_dict = {
        '1': tilt1,
        '5': tilt5,
        }

twist_dict = {
        '1': twist1,
        '5': twist5,
        }

theta = [5, 1]
    
labels = {'1': [r'1$^{\circ}$', '19.6'], 
          '2': [r'2$^{\circ}$', '9.8'], 
          '3': [r'3$^{\circ}$'], 
          '5': [r'5$^{\circ}$', '3.9'],
          '7': [r'7$^{\circ}$', '2.8'], 
          '8': [r'8$^{\circ}$'], 
          '10': [r'10$^{\circ}$', '2.0'], 
          '15': [r'15$^{\circ}$'],}

colors = ['xkcd:dark cyan','xkcd:grass']
i=0
for t in ['5', '1']:
    plt.plot(tilt_dict[t][0] / (tilt.vs * tilt.k_max), tilt_dict[t][1],\
            color = colors[i], linestyle = 'dotted')
    plt.plot(twist_dict[t][0] / (twist.vs * twist.k_max), twist_dict[t][1],\
            color = colors[i])
    i = i+1
    
plt.plot([], [], color = 'xkcd:black', label = 'twist')
plt.plot([], [], color = 'xkcd:black', linestyle = ':', label = 'tilt')

mpl.rcParams['font.size'] = '12'

plt.xlabel('$\omega / \omega_{\mathrm{max}}$')
plt.ylabel(r'$\tau \; \mathrm{(ns)}$')

ax.text(0.13, 0.14, r'$\theta$ = 5$^{\circ}$\\ D = 3.9 nm', verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes, color = 'xkcd:dark cyan', fontsize = 12, weight = 'bold')
ax.text(0.32, 0.77, r'$\theta$ = 1$^{\circ}$\\ D = 19.6 nm', verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes, color = 'xkcd:grass', fontsize = 12, weight = 'bold')

plt.legend()

plt.savefig('tilt_twist_comp.pdf', bbox_inches = 'tight')
    
'''
Log-log Plot
'''
plt.figure()
i=0
for t in ['5', '1']:
    plt.loglog(tilt_dict[t][0] / (tilt.vs * tilt.k_max), tilt_dict[t][1],\
            color = colors[i], label = labels[t][0] + r'; ' + labels[t][1] + ' nm', linestyle = 'dotted')
    plt.loglog(twist_dict[t][0] / (twist.vs * twist.k_max), twist_dict[t][1],\
            color = colors[i], label = labels[t][0] + r'; ' + labels[t][1] + ' nm')
    i = i+1

