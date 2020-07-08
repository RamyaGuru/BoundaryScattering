#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:38:43 2020

@author: ramyagurunathan

Scattering Results Plotting Methods
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ThermalTransport as TT
from ArrayScattering import ArrayScattering as AS
import math
#import AngularVs

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

def diffraction_plot(gb : AS, k_norm_list, Gamma_GBS_list, save = False):
    '''
    plot diffraction peaks
    '''
    plt.figure()
    plt.xlabel(r'$k/k_{\mathrm{max}}$', fontsize=16)
    plt.ylabel(r'$\Gamma \; \mathrm{(ns^{-1})}$', fontsize=16)
    plt.xlim([0,1])
    plt.ylim([0,100])
    plt.plot(k_norm_list, Gamma_GBS_list)
    if save:
        plt.savefig(gb['geom'] + '_' + str(gb['theta']) + '_diffraction.pdf', bbox_inches = 'tight')


        
def convergence_tau_plot(gb : AS, Gamma, n_angle, T, save = False):
    '''
    Converge relaxation time w.r.t the differential angle (d_angle = 2*pi/n_angle)
    used when integrating over incident phonon directions    
    '''
    n_angle_list = np.arange(100, n_angle, 10)
    tau_nlist = []
    for n_angle in n_angle_list:
        tau_nlist.append(TT.tau_spectral(Gamma, gb, gb.k_max / 20., n_angle, T)) 
        print(tau_nlist)
    plt.figure()
    plt.xlabel('n', fontsize=16)
    plt.ylabel(r'$\tau(k) \; \mathrm{(ns)}$', fontsize=16)
    plt.plot(n_angle_list, tau_nlist)
    plt.show(block=False)
    if save == True:
        np.save(str(gb.geom) + str(gb.theta) + 'tau_conv.npy', np.array([n_angle_list, tau_nlist]))


def directional_tau():
    '''
    Plot of tau versus the incident direction of the phonon
    '''
    

def spectral_plots(gb : AS, spectral : dict, prop_list = ['tau', 'transmissivity', 'TBC', 'kappa'], save = False):
    '''
    Plots of spectral properties
    prop_list: specifies properties which should be plotted
    '''
    labels = {
            'tau' : r'$\tau \; \mathrm{(ns)}$',
            'TBC' :  r'$R_K$ $(m^2K/W)$',
            'transmissivity' : 'Transmissivity',
            'kappa' : r'$\kappa_\mathrm{L} \; \mathrm{(W/m/K)}$'
                                        }
    if any(prop_list) not in spectral.keys():
        ValueError()
    for prop in prop_list:
        plt.figure()
        plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
        plt.ylabel(labels[prop], fontsize=16)
        plt.plot([ f / gb.vs for f in spectral['omega']], spectral[prop])
        if save: 
            np.save(str(gb.geom) + str(gb.theta) + 'spectral_update' + prop + '.npy', np.array([spectral['omega'], spectral[prop]]))
            plt.savefig(str(gb.geom) + str(gb.theta) + 'spectral_update' + prop + '.pdf' ,\
                        dpi=400, bbox_inches = 'tight')
        
#def directional_tau_plots(tau_list: list, n_angle):
#    '''
#    Heatmap on a sphere of tau values
#    '''
#    fig = plt.figure()
#    ax = fig.add_subplot( 1, 1, 1, projection='3d')
#    d_angle = (math.pi / 2) / n_angle
#    u = np.arange(d_angle, math.pi / 2 + d_angle, d_angle)
#    v = np.arange(d_angle, math.pi / 2, d_angle)
#    # create the sphere surface
#    XX = 10 * np.outer( np.cos( u ), np.sin( v ) )
#    YY = 10 * np.outer( np.sin( u ), np.sin( v ) )
#    ZZ = 10 * np.outer( np.ones( np.size( u ) ), np.cos( v ) )
#    WW = XX.copy()
#    for i in len(u):
#        WW[i,:] = tau_list[]
#        
        
    

    