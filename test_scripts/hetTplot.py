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

het = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/heterointerfaceNonespectral_updatetau.npy')

vs = (6084 + 5400) / 2

het[1] = het[1]


V = 2E-29
N = 2

k_max = (6 * math.pi**2 / (V * N))**(1 / 3)
dk = k_max / 100

omega = vs * np.arange(dk, k_max, dk)
ax = plt.axes()



'''
Input values for Silicon and Germanium
'''

input_dict = {
        'avg_vs': (6084 + 5400) / 2, #currently using average velocity between material 1 and material 2
             'atmV': [1.97E-29, 2.27E-29],
             'N': 2,
             'bulkmod' : 97.83,
             'shearmod' : 60,
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
                     geom = geom, bvK = True)



'''
Spectral Relaxation Time
'''
het_AMM = np.ones(99) * het[1][0]
#omega_max = (2 / math.pi) * (vs * k_max)
omega_max = vs * k_max

plt.plot(het[0] / hetint.omegaD, het[1], color = 'xkcd:darkish blue', label = 'Full model')
plt.xlabel('$\omega / \omega_{\mathrm{max}}$', fontsize=16)
plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)

#AMM Only
plt.plot(het[0] / hetint.omegaD, het_AMM, color = 'xkcd:darkish blue', linestyle = ':', label = 'Acoustic mismatch only')

mpl.rcParams['font.size'] = '12'
plt.legend(loc = 'bottom left',  frameon = False)

ax = plt.gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.get_minor_formatter().set_scientific(False)
ax.yaxis.get_minor_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
ax.set_ylim(0.01,0.03)
plt.savefig('hetint_tauspectral.pdf', bbox_inches = 'tight')


'''
Thermal Boundary Resistance
'''

tbc = []
wtau = []
Trange = [100, 150, 200, 300]
n_k = 100
dk = hetint.k_max / n_k
k_mags = np.arange(dk, hetint.k_max, dk)

for T in [100, 150, 200, 300]:
    transport = TT.transport_coeffs_from_tau(hetint, k_mags, het[1], T)
    tbc.append((1 / transport['TBC']) * 1E9)
    wtau.append(transport['wtd_tau'])
mpl.rcParams['font.size'] = '16'

plt.figure()
plt.plot(Trange, tbc, color = 'xkcd:darkish blue', label = 'Full model')
plt.xlabel(r'T (K)', fontsize=16)
plt.ylabel(r'$R_K$  (10$^{-9}$ m$^2$K/W)', fontsize=16)

ax = plt.gca()
#ax.set_yscale('log')
#ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.yaxis.get_minor_formatter().set_scientific(False)
#ax.yaxis.get_minor_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))
#plt.savefig('hetint_tbc.pdf', bbox_inches = 'tight')


#ax.yaxis.get_minor_formatter().set_scientific(False)
#plt.yscale('log')

'''
Heterointerface with just the AMM term : TBC
'''

tbc_AMM = []
Trange = [100, 150, 200, 300]

for T in Trange:
    transport = TT.transport_coeffs_from_tau(hetint, k_mags, het_AMM, T)
    tbc_AMM.append((1 / transport['TBC']) * 1E9)
 
print(np.array(tbc_AMM) / np.array(tbc))
'''
Heterointerface with just AMM term : Russian doll tau
'''
wtau_AMM = []
Trange = [100, 150, 200, 300]

for T in Trange:
    transport = TT.transport_coeffs_from_tau(hetint, k_mags, het_AMM, T)
    wtau_AMM.append(transport['wtd_tau'])

print(np.array(wtau) / np.array(wtau_AMM))

plt.plot(Trange, tbc_AMM, linestyle = ':', color = 'xkcd:darkish blue', label = 'Acoustic mismatch only')
plt.xlabel(r'T (K)', fontsize=16)
plt.ylabel(r'$R_K$  (10$^{-9}$ m$^2$K/W)', fontsize=16)
plt.ylim([0,6])
ax = plt.gca()
#ax.set_yscale('log')
#ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.yaxis.get_minor_formatter().set_scientific(False)
#ax.yaxis.get_minor_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f'))
#plt.savefig('hetint_tbc.pdf', bbox_inches = 'tight')

'''
Plot the group velocity
'''
vg = np.zeros(n_k - 1)
i = 0
for k in k_mags:
    vg[i] = hetint.vg_kmag(k)
    i = i+1

plt.figure()
plt.plot(k_mags, vg)