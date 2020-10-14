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
for T in [100, 150,200,250,300]:
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

theta = [10, 7, 3, 1]
colors = []
for t in theta:
    colors.append(scalarMap.to_rgba(t))

fig, ax = plt.subplots(2, sharex = True)
 
i = 0
for tbc in [tbc10, tbc7, tbc3, tbc1]:
    ax[0].plot([100,150,200,250,300], tbc, '-o', color = colors[i], label = labels[i])
    i = i+1

fig.text(0, 0.5, r'$R_K$  (10$^{-9}$ m$^2$K/W)', va='center', rotation='vertical')
ax[1].set_xlabel(r'T (K)')
ax[0].legend(prop = {'size': 12}, loc = (1.02, 0.2))

#Load the experimental data 
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
    ax[1].plot(data[theta][0], data[theta][1],'-o', color = e_colors[j], label = labels[j])
    j = j+1

ax[1].legend(prop = {'size': 12}, loc = (1.02, 0.4))
ax[1].yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[1].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(100))
ax[1].xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
ax[1].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax[1].ticklabel_format(style="plain")


fig.savefig('rk_versusT_pred_exp.pdf', bbox_inches = 'tight')


'''
Spectral Tau Plots
'''

plt.rcParams["figure.figsize"] = [5, 3]

mpl.rcParams['font.size'] = '14'
plt.figure()
ax = plt.gca()


twist_dict = {
        '1': twist1,
        '2': twist2,
        '3': twist3,
        '5': twist5,
        '7': twist7,
        '8' : twist8,
        '10': twist10,
        '15': twist15,
        }
theta = [10, 7, 5, 2, 1]
colors = []
for t in theta:
    colors.append(scalarMap.to_rgba(t))
    
labels = {'1': [r'1$^{\circ}$', '19.6'], 
          '2': [r'2$^{\circ}$', '9.8'], 
          '3': [r'3$^{\circ}$'], 
          '5': [r'5$^{\circ}$', '3.9'],
          '7': [r'7$^{\circ}$', '2.8'], 
          '8': [r'8$^{\circ}$'], 
          '10': [r'10$^{\circ}$', '2.0'], 
          '15': [r'15$^{\circ}$'],}

i=0
for t in ['10', '7', '5', '2', '1']:
    plt.plot(twist_dict[t][0] / (twist.vs * twist.k_max), twist_dict[t][1],\
             color = colors[i], label = labels[t][0] + r'; ' + labels[t][1] + ' nm')
    i = i+1

ax.text(0.85, 0.95, r'$\mathrm{\theta}$; $D$', verticalalignment = 'center', horizontalalignment = 'center',\
        transform = ax.transAxes, color = 'xkcd:black', fontsize = 10)    
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

plt.xlabel('$\omega / \omega_{\mathrm{max}}$')
plt.ylabel(r'$\tau \; \mathrm{(ns)}$')
mpl.rcParams['font.size'] = '10'
#plt.colorbar(scalarMap, shrink = 0.7)
plt.legend(loc = 'upper right', bbox_to_anchor = (1, 0.95), frameon = False)
plt.savefig('twisttauspectral.pdf', bbox_inches = 'tight')



'''
LOGLOG PLOT
'''

mpl.rcParams['font.size'] = '14'
plt.figure()
ax = plt.gca()
i=0
for t in ['10', '7', '5', '2', '1']:
    plt.loglog(twist_dict[t][0] / (twist.vs * twist.k_max), twist_dict[t][1], color = colors[i], label = labels[t])
    i = i+1


'''
Plot the trend lines
'''
plt.loglog(twist10[0][50:] / (twist.vs * twist.k_max), twist10[1][50]*(twist10[0][50:]/twist10[0][50])**(-1.1) * 0.75, color = 'xkcd:black')
#plt.loglog(twist1[0][50:] / (twist.vs * twist.k_max), twist1[1][50]*(twist1[0][50:]/twist1[0][50])**(-1) * 1.4, color = 'xkcd:black')
mpl.rcParams['font.size'] = '10'
#plt.colorbar(scalarMap, shrink = 0.7)
#plt.legend()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(False)
ax.yaxis.get_major_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

plt.xlabel('$\omega / \omega_{\mathrm{max}}$')
plt.ylabel(r'$\tau \; \mathrm{(ns)}$')

plt.savefig('loglog_twisttau.pdf', bbox_inches = 'tight')

'''
Plot the Read Shockley Energy versus the interfacial thermal resistance
'''

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap

trunc = truncate_colormap(plt.get_cmap("viridis_r"), 0, 0.5)
ctNorm  = mcolors.Normalize(vmin=1, vmax=7)
scaMap_trunc = cmx.ScalarMappable(norm=ctNorm, cmap=trunc) 

tbc_list = [tbc1[-1], tbc2[-1], tbc5[-1],  tbc8[-1]]

theta = [1, 2, 5, 7]
colors = []

for t in theta:
    colors.append(scaMap_trunc.to_rgba(t))

GBenergy = []

for t in theta:
    as_obj = AS(**input_dict, geom = 'tilt', theta = t, ax = 1, d_GS = 350E-9, bvK = True)
    print(as_obj.D)
    GBenergy.append(AS.gb_energy(as_obj))

plt.rcParams["figure.figsize"] = [5, 3]
plt.figure()
plt.scatter([tbc for tbc in tbc_list], GBenergy, c = colors, s = 40)
mpl.rcParams['font.size'] = '10'
cbar = plt.colorbar(scaMap_trunc, shrink = 0.5, pad = -0.1, ticks = [1,3,5,7])
#cbar.ax.locator_params(nbins = 4)
cbar.ax.set_yticklabels([r'1$^{\circ}$', r'3$^{\circ}$', r'5$^{\circ}$', r'7$^{\circ}$'])

ax = plt.gca()
#ax.set_aspect(40)
plt.xlabel(r'$R_K$ (m$^2$K/GW)')
plt.ylabel(r'Boundary Energy (J/m$^2)$')
ax.text(0.91, 0.14, r'$\mathrm{\theta}$', verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes, color = 'xkcd:black', fontsize = 20, weight = 'bold')
plt.savefig('TBR_energy.pdf', bbox_inches = 'tight')
    