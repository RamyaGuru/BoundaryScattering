#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 20:41:09 2020

@author: ramyagurunathan

Twist Plot versus Temperature 

Trying the log plot
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
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec


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
Load relaxation times
'''

twist1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist1spectral_updatetau.npy')
twist2 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist2spectral_updatetau.npy')
twist3 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist3spectral_updatetau.npy')
twist5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist5spectral_updatetau.npy')
twist7 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist7spectral_updatetau.npy')
twist8 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist8spectral_updatetau.npy')
twist10 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist10spectral_updatetau.npy')
twist15 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist15spectral_updatetau.npy')

'''
Calculate transport coefficients
'''

labels = [ r'10$^{\circ}$', r'7$^{\circ}$', r'3$^{\circ}$', r'1$^{\circ}$']

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = AS(**input_dict, geom = 'twist', theta = 5, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

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
Set up colormap
'''
viridis = cm = plt.get_cmap('viridis_r') 
cNorm  = mcolors.Normalize(vmin=0, vmax=16)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)   

'''
Load the MD data
'''
fig = plt.figure()


bax =  brokenaxes(xlims=((75,305),(475, 505)), yscale = 'log')

MD = [[500, 500, 500], [0.61, 0.76, 1.1]]

md_color = scalarMap.to_rgba(11.42)

m_list = ['^', 'x', '+']

label_list = ['11.42$^{\circ}$', 'Ref. 6', 'Ref. 21', 'Ref. 4']

bax.scatter([0], [0], marker = 'None', label = '11.42$^{\circ}$')

for m in range(len(MD[0])):
    bax.loglog(MD[0][m], MD[1][m], color = md_color, marker = m_list[m], label = label_list[m+1])
 

theta = [10, 7, 3, 1]
colors = []
for t in theta:
    colors.append(scalarMap.to_rgba(t))

#fig, ax = plt.subplots(2, sharex = True)


i = 0
for tbc in [tbc10, tbc7, tbc3, tbc1]:
    bax.loglog([100,150,200,250,300, 500], tbc, '-o', color = colors[i], label = labels[i])
    i = i+1

fig.text(0, 0.5, r'$R_K$  (10$^{-9}$ m$^2$K/W)', va='center', rotation='vertical')
fig.text(0.5, 0, r'T (K)', va='center')
#ax1.set_xlabel(r'T (K)')
bax.legend(prop = {'size': 12}, loc = (1.02, 0.1))

'''
Load the experimental data
'''
SiSi_twist = np.loadtxt('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/qinghao_sige_rk/Si_Si_twist_datasets.csv', delimiter = ",", skiprows = 2)

angles = ['86.5', '70.4', '15.4', '6.9' ,'3.4']

data = {}

#ax2 = ax1.twinx()

i = 0
for c in range(0,SiSi_twist.shape[1],2):
    data[angles[i]] = [SiSi_twist[:,c], SiSi_twist[:,c+1]]
    i = i+1

e_theta = [6.9, 3.4]
labels = [r'6.9$^{\circ}$', '3.4$^{\circ}$']
e_colors = []
for t in e_theta:
    e_colors.append(scalarMap.to_rgba(t))
  
j = 0
for theta in ['6.9', '3.4']:
    bax.loglog(data[theta][0], data[theta][1],'-o', color = e_colors[j], label = labels[j])
    j = j+1


bax.legend(prop = {'size': 12}, loc = (1.02, 0.4))
fig.savefig('rk_T_pred_exp_MD.pdf', bbox_inches = 'tight')


bax.ticklabel_format(style="plain")




#fig.savefig('rk_T_pred_exp_MD_log.pdf', bbox_inches = 'tight')

