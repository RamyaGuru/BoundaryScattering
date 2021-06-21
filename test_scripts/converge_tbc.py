#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 23:18:10 2021

@author: ramyagurunathan

Converge the thermal boundary resistance
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

'''
Load relaxation times

'''

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'shearmod' : 60E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = AS(**input_dict, geom = 'twist', theta = 11, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)


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


#def drop_some_values(tau_list, n):
#    '''
#    Drops every n element in tau_list
#    '''
#    return tau_list[::n]



def convergence_plot_wrt_n(tau_list, del_steps : int, n : int):
    tau_dict = {}
    tau_dict['0'] = tau_list
    red_tau = tau_list
    for d in range(del_steps):
        red_tau = [red_tau[0][1::n], red_tau[1][1::n]]
        tau_dict[str(d + 1)] = red_tau
    return tau_dict

T = 300
x_list = []
y_list = []

tau_dict = convergence_plot_wrt_n(twist1, del_steps = 4, n = 2)
for k in tau_dict.keys():
    x_list.append(len(tau_dict[k][0]))
    n_k = len(tau_dict[k][0])
    dk = twist.k_max / n_k
    k_mags = np.arange(dk, twist.k_max, dk)
    transport = TT.transport_coeffs_from_tau(twist, k_mags, tau_dict[k][1], T)
    y_list.append(transport['TBC'])

plt.scatter(x_list, y_list)
    
