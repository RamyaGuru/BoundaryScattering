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
    plt.plot(k_norm_list, Gamma_GBS_list)
    if save:
        plt.savefig(gb['geom'] + '_' + str(gb['theta']) + '_diffraction.pdf', bbox_inches = 'tight')


        
def convergence_tau_plot(gb : AS, Gamma, n_angle, T):
    '''
    Converge relaxation time w.r.t the differential angle (d_angle = 2*pi/n_angle)
    used when integrating over incident phonon directions    
    '''
    n_angle_list = np.arange(4, n_angle, 2)
    tau_nlist = []
    for n_angle in n_angle_list:
        tau_nlist.append(TT.tau_spectral(Gamma, gb, gb.k_max / 5., n_angle, T))    
    plt.figure()
    plt.xlabel('n', fontsize=16)
    plt.ylabel(r'$\tau(k)^{-1} \; \mathrm{(ns^{-1})}$', fontsize=16)
    plt.plot(n_angle_list, tau_nlist)
    plt.show(block=False)


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
            plt.savefig(gb['geom'] + gb['theta'] + 'spectral_' + prop ,\
                        dpi=400, bbox_inches = 'tight')
        
    
