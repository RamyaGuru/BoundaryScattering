#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:23:25 2020

@author: ramyagurunathan

Plot the twist boundary versus angle thermal boundary conductance
"""
import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/scattering_scripts')

from ArrayScattering import ArrayScattering
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

#twist = np.load('../datafiles/3_twist_vsT.npz')
#twist1 = np.load('../datafiles/15.4_twist_vsT.npz')
#twist2 = np.load('../datafiles/86.5_twist_vsT.npz')
#twist3 = np.load('70.4_twist_vsT.npz')

color = ['xkcd:light indigo', 'xkcd:tree green', 'xkcd:blood red']
labels = [r'3.4$^{\circ}$',  r'6.9$^{\circ}$', r'15.4$^{\circ}$']

i = 0
fig, ax = plt.subplots(2, sharex = True)

#for tb in [twist, twist1, twist2]:
#    ax[0].plot(twist['arr_0'], tb['arr_2'] * 1e9, 'o-', color = color[i], label = labels[i])
#    i = i+1
    
#ax[0].set_ylabel(r'$R_K$  (10$^{-9}$ m$^2$K/W)'
fig.text(0, 0.5, r'$R_K$  (10$^{-9}$ m$^2$K/W)', va='center', rotation='vertical')
ax[1].set_xlabel(r'T (K)')
ax[0].legend(loc = 'lower left', prop = {'size': 12})
ax[0].yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax[0].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[0].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[0].ticklabel_format(style="plain")

    
'''
New plot with the experimental data 
'''

#Load the experimental data 
SiSi_twist = np.loadtxt('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/qinghao_sige_rk/Si_Si_twist_datasets.csv', delimiter = ",", skiprows = 2)

angles = ['86.5', '70.4', '15.4', '6.9' ,'3.4']

data = {}

#ax2 = ax1.twinx()

i = 0
for c in range(0,SiSi_twist.shape[1],2):
    data[angles[i]] = [SiSi_twist[:,c], SiSi_twist[:,c+1]]
    i = i+1
  
j = 0
for theta in ['15.4', '6.9', '3.4']:
    ax[1].plot(data[theta][0], data[theta][1],'-o', color = color[j], label = labels[2-j])
    j = j+1

ax[1].legend(prop = {'size': 12}, loc = (1.02, 0.4))
ax[1].yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[1].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(100))
ax[1].xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
ax[1].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[1].ticklabel_format(style="plain")
'''
Add the MD simulation data
'''
#Load the MD simulation data
md_11deg = np.genfromtxt('../datafiles/11deg_MD.csv', delimiter = ',')
md_16deg = np.genfromtxt('../datafiles/16deg_MD.csv', delimiter = ',')
md_36deg = np.genfromtxt('../datafiles/36deg_MD.csv', delimiter = ',')
#ax[0].scatter(md_11deg[0,0], md_11deg[0,1], s = 60, color = 'xkcd:plum', marker = 'x', label = r'11.42$^{\circ}$ MD')
#ax[0].scatter(300, 0.61, s = 60, color = 'xkcd:plum', marker = 'x') #Schelling point at 500 K
#ax1.scatter(md_36deg[0,0], md_36deg[0,1], s = 60, color = 'xkcd:black', marker = 'x')

fig.savefig('rk_versusT_pred_exp.pdf', bbox_inches = 'tight')

'''
Add Twist Interface Predictions
'''
input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = ArrayScattering(**input_dict, geom = 'twist', theta = 2, ax = {'n' : 1, 'm' : 2}, d_GS = 350E-9)

'''
Load Relaxation Times
'''
#300 K
twist2 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist2spectral_tau.npy')
twist5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist5spectral_tau.npy')
twist10 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist10spectral_tau.npy')
twist15 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist15spectral_tau.npy')

#'''
#Twist: T = 300 K 
#'''
#T = 300
#
#transport2 = TT.transport_coeffs_from_tau(twist, twist2[0] / twist.vs, twist2[1], T)
#transport5 = TT.transport_coeffs_from_tau(twist, twist5[0] / twist.vs, twist5[1], T)
#transport10 = TT.transport_coeffs_from_tau(twist, twist10[0] / twist.vs, twist10[1], T)
#transport15 = TT.transport_coeffs_from_tau(twist, twist15[0] / twist.vs, twist15[1], T)
#
color = ['xkcd:black', 'xkcd:plum', 'xkcd:blood red', 'xkcd:dark salmon']
labels = [r'15$^{\circ}$', r'10$^{\circ}$', r'5$^{\circ}$', r'2$^{\circ}$']
##i = 0
##for tbc in [transport15['TBC'], transport10['TBC'], transport5['TBC'], transport2['TBC']]:
##    ax[0].scatter(300, (1/tbc)*1e9, color = color[i], label = labels[i])
##    i = i+1
##
##ax[0].legend(prop = {'size': 12})
##ax[0].scatter([300,300,300,300], [(1/tbc)*1e9 for tbc in [transport2['TBC'], transport5['TBC'], transport10['TBC'], transport15['TBC']]])
#
#'''
#Twist: T = 200 K
#'''
#T = 200
#twist15_200 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist15spectral_tauT200.npy')
#transport15_200 = TT.transport_coeffs_from_tau(twist, twist15_200[0] / twist.vs, twist15_200[1], 200)
#
##ax[0].scatter(200, (1/transport15_200['TBC'])*1E9, color = 'xkcd:black')
#
#'''
#Twist: T = 150 K
#'''
#T = 150
#twist15_150 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist15spectral_tauT150.npy')
#transport15_150 = TT.transport_coeffs_from_tau(twist, twist15_150[0] / twist.vs, twist15_150[1], 150)
#       
#tbc = {'15': [[150, 200, 300], [transport15_150['TBC'], transport15_200['TBC'], transport15['TBC']]],\
#       '10' : [[300], [transport10['TBC']]],\
#       '5' : [[300],[transport5['TBC']]],\
#       '2' : [[300],[transport2['TBC']]]}

'''
Calculate Transport Properties at a Range of Temperatures
'''
tbc2 = []
tbc5 = []
tbc10 = []
tbc15 = []

for T in [100, 150, 200, 300]:
    transport2 = TT.transport_coeffs_from_tau(twist, twist2[0] / twist.vs, twist2[1], T)
    transport5 = TT.transport_coeffs_from_tau(twist, twist5[0] / twist.vs, twist5[1], T)
    transport10 = TT.transport_coeffs_from_tau(twist, twist10[0] / twist.vs, twist10[1], T)
    transport15 = TT.transport_coeffs_from_tau(twist, twist15[0] / twist.vs, twist15[1], T)
    tbc2.append((1 / transport2['TBC']) * 1E9)
    tbc5.append((1 / transport5['TBC']) * 1E9)
    tbc10.append((1 / transport10['TBC']) * 1E9)
    tbc15.append((1 / transport15['TBC']) * 1E9)
    
i = 0
for tbc in [tbc15, tbc10, tbc5, tbc2]:
    ax[0].plot([100, 150, 200, 300], tbc, '-o', color = color[i], label = labels[i])
    i = i+1

ax[0].legend(prop = {'size': 12}, loc = (1.02,0.2))
fig.savefig('TBC_comp.pdf', bbox_inches = 'tight')