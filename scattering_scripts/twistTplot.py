#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:23:25 2020

@author: ramyagurunathan

Plot the twist boundary versus angle thermal boundary conductance
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcdefaults()
mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = '18'
mpl.rcParams['text.usetex'] = True

mpl.rcParams['mathtext.fontset'] = 'custom'

mpl.rcParams['mathtext.bf'] = 'Apple Symbols'

twist = np.load('../datafiles/3_twist_vsT.npz')
twist1 = np.load('../datafiles/15.4_twist_vsT.npz')
twist2 = np.load('../datafiles/86.5_twist_vsT.npz')
#twist3 = np.load('70.4_twist_vsT.npz')

color = ['xkcd:light indigo', 'xkcd:tree green', 'xkcd:blood red']
labels = [r'3$^{\circ}$', r'15.4$^{\circ}$', r'86.4$^{\circ}$']

i = 1
fig, ax1 = plt.subplots()
for tb in [twist1, twist2]:
    ax1.loglog(twist['arr_0'], tb['arr_2'] * 1e9, 'o', color = color[i], label = labels[i])
    i = i+1
    
ax1.set_ylabel(r'$R_K$  (10$^{-9}$ m$^2$K/W)')
ax1.set_xlabel(r'T (K)')
ax1.legend()


#
#j=0
#plt.figure()
#for tc in [twist, twist1, twist2]:
#    plt.plot(twist['arr_0'], tc['arr_1'], color = color[j])
#    j = j+1
    
'''
New plot with the experimental data 
'''
plt.figure()
#Load the experimental data 
SiSi_twist = np.loadtxt('../../../qinghao_sige_rk/Si_Si_twist_datasets.csv', delimiter = ",", skiprows = 2)

angles = ['86.5', '70.4', '15.4', '6.9' ,'3.4']

data = {}

ax2 = ax1.twinx()

i = 0
for c in range(0,SiSi_twist.shape[1],2):
    data[angles[i]] = [SiSi_twist[:,c], SiSi_twist[:,c+1]]
    i = i+1
  
j = 1
for theta in ['15.4','86.5']:
    ax2.loglog(data[theta][0], data[theta][1], linestyle = ':', color = color[j])
    j = j+1


'''
Add the MD simulation data
'''
#Load the MD simulation data
md_11deg = np.genfromtxt('../datafiles/11deg_MD.csv', delimiter = ',')
md_16deg = np.genfromtxt('../datafiles/16deg_MD.csv', delimiter = ',')
md_36deg = np.genfromtxt('../datafiles/36deg_MD.csv', delimiter = ',')
#ax1.scatter(md_16deg[0,0], md_16deg[0,1], s = 60, color = 'xkcd:black', marker = 'x')
#ax1.scatter(md_36deg[0,0], md_36deg[0,1], s = 60, color = 'xkcd:black', marker = 'x')

plt.savefig('rk_versusT_pred_exp.pdf', bbox_inches = 'tight')