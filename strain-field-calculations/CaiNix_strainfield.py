#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 20:32:54 2020

@author: ramyagurunathan

Peierls Nabarro

Heterointerface Strain field

From: Cai and Nix Textbook
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, sinh, cosh
from math import pi
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib as mpl

mpl.rcParams['font.size'] = '14'
mpl.rcParams.update({'figure.autolayout': True})
'''
Inputs
'''
a = 5.4e-10
a2 = 5.658e-10 
V = 2E-29
N = 2

bx = 0
by = .34#(V * N) ** (1 / 3) 

P = a/(a2-a)
p = 2.72#P * b

mu = 6.81e10 
nu = 0.27

def X(x):
    return x / p

def Y(y):
    return y / p

def sigma_xx(x, y):
    term1 = -bx * sin(2 * pi * Y(y)) * ((2 * pi * X(x) * sinh(2 * pi * X(x)) +\
                     cosh(2 * pi * X(x))- cos(2 * pi * Y(y))) /\
                      (cosh(2 * pi * X(x)) - cos(2 * pi * Y(y)))**2)
    term2 = by * 2 * pi * X(x) * ((cosh(2 * pi * X(x)) * cos(2 * pi * Y(y)) - 1) /\
                           (cosh(2 * pi * X(x)) - cos(2 * pi * Y(y)))**2)
    return (mu / (2 * (1 - nu) * p)) * (term1 + term2)

def sigma_yy(x, y):
    term1 = bx * sin(2 * pi * Y(y)) * (2 * pi * X(x) * sinh(2 * pi * X(x)) -\
                     cosh(2 * pi * X(x)) + cos(2 * pi * Y(y))) /\
                     (cosh(2 * pi * X(x)) - cos(2 * pi * Y(y)))**2 
    term2 = -by * (2 * pi * X(x) * (cosh(2 * pi * X(x)) * cos(2 * pi * Y(y)) - 1) -\
                 2 * sinh(2 * pi * X(x)) * (cosh(2 * pi * X(x)) - cos(2 * pi * Y(y))))/\
                           (cosh(2 * pi * X(x)) - cos(2 * pi * Y(y)))**2  
    return (mu / (2 * (1 - nu) * p)) * (term1 + term2)

def sigma_yx(x, y):
    term1 = bx * 2 * pi * X(x) * ((cosh(2 * pi * X(x)) * cos(2 * pi * Y(y)) - 1) /\
                           (cosh(2 * pi * X(x)) * cos(2 * pi * Y(y)))**2)
    term2 = by * sin(2 * pi * Y(y)) * (2 * pi * X(x) * sinh(2 * pi * X(x)) +\
                     cosh(2 * pi * X(x))- cos(2 * pi * Y(y)))
    return (mu / (2 * (1 - nu) * p)) * (term1 + term2)

def eps_yy(x,y):
    sigyy = sigma_yy(x,y)
    sigxx = sigma_xx(x,y)
    epsyy = (1 / (2 * mu))*((1 - nu) * sigyy - nu * sigxx)
    return epsyy

if __name__ == '__main__':
#    '''
#    Stress Plots
#    '''
#    xlim = 6
#    ylim = 6
#    nelem = 1000
#    x = np.linspace(-xlim,xlim,nelem)
#    y = np.linspace(-ylim,ylim,nelem)
#    Xax, Yax = np.meshgrid(x,y)
#    sigmaxx = np.zeros([0, nelem])
#    for why in y:
#        sigmaxx_zet = sigma_yy(x, why)
#        sigmaxx_zet = np.transpose(sigmaxx_zet.reshape(nelem,1))
#        sigmaxx = np.append(sigmaxx, sigmaxx_zet, axis = 0)
#    
#    fig = plt.figure()
#    ax = Axes3D(fig)
#    ax = plt.axes(projection='3d')
#    sigmaxx[sigmaxx > 0.5 * mu] = 0.5 * mu
#    sigmaxx[sigmaxx < -0.5 * mu] = -0.5 * mu
#    ax.plot_surface(Yax, Xax, sigmaxx / mu, cmap = 'coolwarm')
#    ax.set_zlim3d(-0.5, 0.5)
#    ax.view_init(elev = 20, azim = -20)
#    ax.set_xlabel('y (nm)')
#    ax.set_ylabel('x (nm)', labelpad = 10)
#    ax.set_zlabel(r'$\sigma_{yy}/\mu$')
#    plt.tight_layout()
#    plt.savefig('misfit_disl_schem.pdf', bbox_inches = 'tight')
#    '''
#    Cross-section figure
#    '''
#    plt.figure()
#    plt.plot(x, sigmaxx[120,:] / mu, 'xkcd:black')
#    plt.xlabel('x')
#    plt.ylabel(r'$\sigma_{yy}/\mu$')
#    plt.savefig('misfit_disl_xsection.pdf', bbox_inches = 'tight')
    
    '''
    Strain Plots
    '''
    xlim = 4
    ylim = 6
    nelem = 1000
    x = np.linspace(-xlim,xlim,nelem)
    y = np.linspace(-ylim,ylim,nelem)
    Xax, Yax = np.meshgrid(x,y)
    epsxx = np.zeros([0, nelem])
    for why in y:
        epsxx_zet = eps_yy(x, why)
        epsxx_zet = np.transpose(epsxx_zet.reshape(nelem,1))
        epsxx = np.append(epsxx, epsxx_zet, axis = 0)
    
    fig = plt.figure()
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = [4.5, 3]
    ax = Axes3D(fig)
    ax = plt.axes(projection='3d')
    epsxx[epsxx > 0.2] = 0.2
    epsxx[epsxx < -0.2] = -0.2
    ax.plot_surface(Yax, Xax, epsxx, cmap = 'coolwarm')
    ax.set_zlim3d(-0.2, 0.2)
    ax.view_init(elev = 20, azim = -30)
    ax.set_xlabel('y (nm)')
    ax.set_ylabel('x (nm)', labelpad = 10)
    ax.set_zlabel(r'$\epsilon_{yy}$')
    plt.savefig('misfit_disl_strain.pdf', bbox_inches = 'tight')
    '''
    Cross-section figure
    '''
    plt.figure()
    plt.plot(x, epsxx[120,:], 'xkcd:black')
    plt.xlabel('x')
    plt.ylabel(r'$\epsilon_{yy}$')
    plt.savefig('misfit_disl_xc_strain.pdf', bbox_inches = 'tight')    