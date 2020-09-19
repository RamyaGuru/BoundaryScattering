#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 09:53:26 2020

@author: ramyagurunathan

Heterointerface Plots
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/scattering_scripts')

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import HetIntScattering as HS

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import math

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
plt.rcParams["figure.figsize"] = [5, 3]

mpl.rcParams['mathtext.fontset'] = 'custom'

mpl.rcParams['mathtext.bf'] = 'Apple Symbols'

het = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020/hetintspectral_updatetau.npy')


vs = (6084 + 5400) / 2
V = 2E-29
N = 2

k_max = (6 * math.pi**2 / (V * N))**(1 / 3)
dk = k_max / 100

omega = vs * np.arange(dk, k_max, dk)
ax = plt.axes()

'''
Spectral Relaxation Time
'''

plt.plot(het[0] / (vs * k_max), het[1], color = 'xkcd:darkish blue')
plt.xlabel('$\omega / \omega_{\mathrm{max}}$', fontsize=16)
plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)

ax = plt.gca()
#ax.set_xscale('log')
ax.set_yscale('log')
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.get_minor_formatter().set_scientific(False)
ax.yaxis.get_minor_formatter().set_useOffset(False)
ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
plt.savefig('hetint_tauspectral.pdf', bbox_inches = 'tight')

'''
Thermal Boundary Resistance
'''

'''
Input values for Silicon and Germanium
'''

input_dict = {
        'avg_vs': (6084 + 5400) / 2, #currently using average velocity between material 1 and material 2
             'atmV': [1.97E-29, 2.27E-29],
             'N': 2,
             'bulkmod' : 97.83,
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

hetint = HS.initialize(input_dict, cmat = [cmat1, cmat2], density = [density1, density2],\
                     geom = geom)

tbc = []
Trange = [100, 150, 200, 300]

for T in [100, 150, 200, 300]:
    transport = TT.transport_coeffs_from_tau(hetint, het[0] / hetint.vs, het[1], T)
    tbc.append((1 / transport['TBC']) * 1E9)

plt.figure()
plt.plot(Trange, tbc, color = 'xkcd:darkish blue')
plt.xlabel(r'T (K)', fontsize=16)
plt.ylabel(r'$R_K$  (10$^{-9}$ m$^2$K/W)', fontsize=16)

ax = plt.gca()
#ax.set_yscale('log')
#ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.yaxis.get_minor_formatter().set_scientific(False)
#ax.yaxis.get_minor_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))
plt.savefig('hetint_tbc.pdf', bbox_inches = 'tight')


#ax.yaxis.get_minor_formatter().set_scientific(False)
#plt.yscale('log')

