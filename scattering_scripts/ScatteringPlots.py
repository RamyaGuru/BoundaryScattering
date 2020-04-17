#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:38:43 2020

@author: ramyagurunathan

Scattering Results Plotting Methods
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import ArrayScattering

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

def diffraction_plot(gb : ArrayScattering, k_norm_list, Gamma_GBS_list, save = False):
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
    

#def tau_spectral_plot():
#    '''
#    plot the spectral relaxation time
#    '''
#    
#    
#def kappa_plot():
#    '''
#    plot the spectral thermal conductivity
#    '''
#    