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

twist1 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_update_bvK/twist1spectral_update2tau.npy')
twist2 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist2spectral_updatetau.npy')
twist3 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/lit_angles/twist3.4spectral_update2tau.npy')
twist5 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist5spectral_updatetau.npy')
twist7 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/lit_angles/twist6.9spectral_update2tau.npy')
twist8 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist8spectral_updatetau.npy')
twist10 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_update_bvK/twist10spectral_update2tau.npy')
twist11 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/lit_angles/twist11.42spectral_update2tau.npy')
twist12 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_update_bvK/twist12spectral_update2tau.npy')
twist15 = np.load('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/fall2020_2/twist15spectral_updatetau.npy')

'''
Calculate transport coefficients
'''

labels = [  r'11.42$^{\circ}$', r'6.9$^{\circ}$', r'3.4$^{\circ}$', r'1$^{\circ}$']

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'shearmod' : 60E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = AS(**input_dict, geom = 'twist', theta = 11, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)

tbc1 = []
tbc2 = []
tbc3 = []
tbc5 = []
tbc7 = []
tbc8 = []
tbc10 = []
tbc11 = []
tbc12 = []
tbc15 = []

n_k = 100
dk = twist.k_max / n_k
k_mags = np.arange(dk, twist.k_max, dk)
for T in [100, 150,200,250,300, 500]:
    transport1 = TT.transport_coeffs_from_tau(twist, k_mags, twist1[1], T)
    transport2 = TT.transport_coeffs_from_tau(twist, k_mags, twist2[1], T)
    transport3 = TT.transport_coeffs_from_tau(twist, k_mags, twist3[1], T)
    transport5 = TT.transport_coeffs_from_tau(twist, k_mags, twist5[1], T)
    transport7 = TT.transport_coeffs_from_tau(twist, k_mags, twist7[1], T)
    transport8 = TT.transport_coeffs_from_tau(twist, k_mags, twist8[1], T)
    transport10 = TT.transport_coeffs_from_tau(twist, k_mags, twist10[1], T)
    transport11 = TT.transport_coeffs_from_tau(twist, k_mags, twist11[1], T)
    transport12 = TT.transport_coeffs_from_tau(twist, k_mags, twist12[1], T)
    transport15 = TT.transport_coeffs_from_tau(twist, k_mags, twist15[1], T)
    tbc1.append((1 / transport1['TBC']) * 1E9)
    tbc2.append((1 / transport2['TBC']) * 1E9)
    tbc3.append((1 / transport3['TBC']) * 1E9)
    tbc5.append((1 / transport5['TBC']) * 1E9)
    tbc7.append((1 / transport7['TBC']) * 1E9)
    tbc8.append((1 / transport8['TBC']) * 1E9)
    tbc10.append((1 / transport10['TBC']) * 1E9)
    tbc11.append((1 / transport11['TBC']) * 1E9)
    tbc12.append((1 / transport12['TBC']) * 1E9)
    tbc15.append((1 / transport15['TBC']) * 1E9)
    
'''
Set up colormap
'''
viridis = cm = plt.get_cmap('viridis_r') 
cNorm  = mcolors.Normalize(vmin=0, vmax=16)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)   

#'''
#Load the MD data
#'''
fig = plt.figure()
sps1, sps2 = GridSpec(2,1)

ax0 =  brokenaxes(xlims=((75,305),(475, 505)), subplot_spec=sps1)
ax1 =  brokenaxes(xlims=((75,305),(475, 505)), subplot_spec=sps2)


MD = [[500, 500, 500], [0.61, 0.76, 1.1]]

md_color = scalarMap.to_rgba(11.42)

m_list = ['^', 'x', '+']

label_list = ['11.42$^{\circ}$', 'Ref. 2', 'Ref. 3', 'Ref. 4']

ax0.scatter([0], [0], marker = 'None', label = '11.42$^{\circ}$') 

for m in range(len(MD[0])):
    ax0.scatter(MD[0][m], MD[1][m], color = md_color, marker = m_list[m], label = label_list[m+1])

leg1 = ax0.legend(prop = {'size': 12}, loc = (1.02, 0.1))

theta = [11, 7, 3, 1]
colors = []
for t in theta:
    colors.append(scalarMap.to_rgba(t))

#fig, ax = plt.subplots(2, sharex = True)


i = 0
for tbc in [tbc11, tbc7, tbc3, tbc1]:
    ax0.semilogy([100,150,200,250,300, 500], tbc, '-o', color = colors[i], label = labels[i])
    i = i+1

fig.text(-0.04, 0.5, r'$R_K$  (10$^{-9}$ m$^2$K/W)', va='center', rotation='vertical')
fig.text(0.5, 0, r'T (K)', va='center')
#ax1.set_xlabel(r'T (K)')
leg2 = ax0.legend(prop = {'size': 12}, loc = (1.02, 0.01))

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
    
#ax1.scatter([0], [0], marker = 'None', label = 'Ref. 5') 
  
j = 0
for theta in ['6.9', '3.4']:
    ax1.plot(data[theta][0], data[theta][1],'-o', color = e_colors[j], label = labels[j])
    j = j+1


ax1.legend(prop = {'size': 12}, loc = (1.02, 0.4))
fig.savefig('rk_T_pred_exp_MD.pdf', bbox_inches = 'tight')


#ax1.ticklabel_format(style="plain")




#fig.savefig('rk_T_pred_exp_MD_bvK.pdf', bbox_inches = 'tight')

#%%

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
        '11':twist11,
        '12':twist12,
        '15': twist15,
        }
theta = [10, 7, 1]
colors = []
for t in theta:
    colors.append(scalarMap.to_rgba(t))
    
labels = {'1': [r'1$^{\circ}$', '19.6'], 
          '2': [r'2$^{\circ}$', '9.8'], 
          '3': [r'3$^{\circ}$', '5'], 
          '5': [r'5$^{\circ}$', '3.9'],
          '7': [r'7$^{\circ}$', '2.8'], 
          '8': [r'8$^{\circ}$'], 
          '10': [r'10$^{\circ}$', '2.0'],
          '11': [r'11$^{\circ}$', '1.8'],
          '12' : [r'12$^{\circ}$'],
          '15': [r'15$^{\circ}$'],}

i=0
for t in ['10', '7', '1']:
    plt.plot(twist_dict[t][0] / (twist.omegaD), twist_dict[t][1],\
             color = colors[i], label = labels[t][0] + r'; ' + labels[t][1] + ' nm')
    i = i+1

ax.text(0.85, 0.75, r'$\mathrm{\theta}$; $D$', verticalalignment = 'center', horizontalalignment = 'center',\
        transform = ax.transAxes, color = 'xkcd:black', fontsize = 10)    
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

plt.xlabel('$\omega / \omega_{\mathrm{max}}$')
plt.ylabel(r'$\tau \; \mathrm{(ns)}$')
mpl.rcParams['font.size'] = '10'
#plt.colorbar(scalarMap, shrink = 0.7)
plt.legend(loc = 'upper right', bbox_to_anchor = (1, 0.75), frameon = False)
#plt.savefig('twisttauspectral.pdf', bbox_inches = 'tight')



'''
LOGLOG PLOT
'''

mpl.rcParams['font.size'] = '14'
plt.figure()
ax = plt.gca()
i=0
for t in ['10', '7', '1']:
    plt.loglog(twist_dict[t][0][:80] / (twist.omegaD), twist_dict[t][1][:80], color = colors[i], label = labels[t])
    i = i+1


'''
Plot the trend lines
'''
plt.loglog(twist10[0][50:] / (twist.omegaD), twist10[1][50]*(twist10[0][50:]/twist10[0][50])**(-1.1) * 0.75, color = 'xkcd:black')
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
plt.ylabel(r'$\tau \; \mathrm{(ns)}$', labelpad = 0)

#plt.savefig('loglog_twisttau.pdf', bbox_inches = 'tight')

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

trunc = truncate_colormap(plt.get_cmap("viridis_r"), 0, 0.6875)
ctNorm  = mcolors.Normalize(vmin=1, vmax=11)
scaMap_trunc = cmx.ScalarMappable(norm=ctNorm, cmap=trunc) 

tbc_list = [tbc1[-1], tbc2[-1],  tbc5[-1],   tbc8[-1], tbc11[-1]]#, tbc10[-1], tbc15[-1]]

theta = [1, 2,  5,  8, 11]# 10, 15]
colors = []

for t in theta:
    colors.append(scalarMap.to_rgba(t))

Strain_energy = []
Strain_energy2 = []
Tot_energy = []
Tot_energy2 = []

for t in theta:
    as_obj = AS(**input_dict, geom = 'tilt', theta = t, ax = 1, d_GS = 350E-9)
    Strain_energy.append(AS.gb_energy(as_obj))
    Tot_energy.append(AS.gb_energy(as_obj))



plt.rcParams["figure.figsize"] = [5, 3]
plt.figure()
plt.scatter([tbc for tbc in tbc_list], Strain_energy, c = colors, s = 40)
mpl.rcParams['font.size'] = '10'
cbar = plt.colorbar(scaMap_trunc, shrink = 0.5, pad = -0.1)
#cbar.ax.locator_params(nbins = 4)
cbar.ax.set_yticklabels([r'1$^{\circ}$', r'3$^{\circ}$', r'5$^{\circ}$', r'7$^{\circ}$',  r'10$^{\circ}$', r'15$^{\circ}$'])

ax = plt.gca()
#ax.set_aspect(40)
plt.xlabel(r'$R_K$ (m$^2$K/GW)')
plt.ylabel(r'Interfacial Energy (J/m$^2)$')
ax.text(0.91, 0.14, r'$\mathrm{\theta}$', verticalalignment = 'center', horizontalalignment = 'center', transform = ax.transAxes, color = 'xkcd:black', fontsize = 20, weight = 'bold')
#plt.savefig('TBR_energy.pdf', bbox_inches = 'tight')


'''
Try Double y-axis plot
'''

fig, ax1 = plt.subplots()

color = 'xkcd:blood red'
ax1.set_xlabel(r'Twist Angle $\theta$ (degrees)')
ax1.set_ylabel('Interfacial Energy (J/m$^2$)')
ax1.plot(theta, Tot_energy, color = color)
#ax1.plot(angles, Tot_energy2, color = color, linestyle = ':')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'xkcd:medium blue'
ax2.set_xlabel(r'Twist Angle $\theta$ (degrees)')
ax2.set_ylabel(r'$R_K$ (m$^2$K/GW)', labelpad = 7)
ax2.scatter(theta, [tbc for tbc in tbc_list], color = color, s = 40)
ax2.tick_params(axis='y', labelcolor=color)

fig.savefig('TBR_energy_bvK.pdf', bbox_inches = 'tight')

'''
Get the AMM contribution to the thermal boundary resistance
'''
twist1_amm = np.ones(100) * twist1[1][0]
twist3_amm = np.ones(100) * twist3[1][0]
twist7_amm = np.ones(100) * twist7[1][0]
twist10_amm = np.ones(100) * twist10[1][0]

for t, tbc in zip([twist1_amm, twist3_amm, twist7_amm, twist10_amm], [tbc1, tbc3, tbc7, tbc10]):
    tbc_AMM = []
    for T in [100, 150,200,250,300, 500]:
        transport = TT.transport_coeffs_from_tau(twist, twist1[0] / twist.vs, t, T)
        tbc_AMM.append((1 / transport['TBC']) * 1E9)
    print(np.array(tbc_AMM) / np.array(tbc))
    
'''
Get the AMM contribution to the Russian Doll relaxation time
'''
twist1_amm = np.ones(100) * twist1[1][0]
twist3_amm = np.ones(100) * twist3[1][0]
twist7_amm = np.ones(100) * twist7[1][0]
twist10_amm = np.ones(100) * twist10[1][0]

for t, wtau in zip([twist1_amm, twist3_amm, twist7_amm, twist10_amm], [trans['wtd_tau'] for trans in [transport1,\
                   transport3, transport7, transport10]]):
    wtau_AMM = []
    for T in [100, 150,200,250,300, 500]:
        transport = TT.transport_coeffs_from_tau(twist, twist1[0] / twist.vs, t, T)
        wtau_AMM.append(transport['wtd_tau'])
    print(np.array(wtau) / np.array(wtau_AMM))
    