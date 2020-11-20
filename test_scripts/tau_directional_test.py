#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 13:50:03 2020

@author: ramyagurunathan

Tau directional Plots
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import math
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

fname = '/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/twist5oldtau_directional_2_hf.pkl'

with open(fname, 'rb') as file:
    tau = pickle.load(file)

tau = np.array(tau)

n_angle = (len(tau))**(1/2)

        
def directional_tau_plots(tau_list: list, n_angle, save = None):
    '''
    Heatmap on a sphere of tau values
    '''
    fig = plt.figure()
    ax = fig.add_subplot( 1, 1, 1, projection='3d')
    d_angle = (math.pi) / n_angle
    u = np.arange(d_angle, math.pi + d_angle, d_angle)
    v = np.arange(d_angle, math.pi, d_angle)
    # create the sphere surface
    XX = 10 * np.outer( np.cos( u ), np.sin( v ) )
    YY = 10 * np.outer( np.sin( u ), np.sin( v ) )
    ZZ = 10 * np.outer( np.ones( np.size( u ) ), np.cos( v ) )
    WW = XX.copy()
    for i in range(len(u)):
        start = i * n_angle
        stop = i * n_angle + (n_angle - 1)
        WW[i,:] = 1/tau_list[start:stop, 2]
#    print(np.amin(WW))
#    WW = WW / np.amax(WW)
    norm = norm = cm.colors.Normalize(vmax=2, vmin=0)
    scamap = plt.cm.ScalarMappable(norm = norm, cmap = 'plasma')
    fcolors = scamap.to_rgba(WW)
    ax.plot_surface(XX, YY, ZZ, cstride = 1, rstride = 1, facecolors = fcolors, cmap = 'plasma')
    cbar = fig.colorbar(scamap, ax = ax, pad = -0.04, shrink = 0.5)
    cbar.set_label(r'$\tau^{-1} (GHz)$', size = 14)
    ax.view_init(azim = 80, elev = 10)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.grid(False)
    plt.tight_layout()
    if save:
        plt.savefig(save + '.pdf', bbox_inches = 'tight')
    
directional_tau_plots(tau, int((len(tau))**(1/2)))
