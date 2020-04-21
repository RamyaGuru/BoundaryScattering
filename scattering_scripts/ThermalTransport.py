#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:37:01 2020

@author: ramyagurunathan

Transport Properties Module
"""

import numpy as np
from ArrayScattering import ArrayScattering as AS
import math

'''
Constants
'''
hbar = 6.6260696 * 10**-34 / (2 * math.pi)
kB = 1.38064852e-23


'''
PROPERTIES: options for calculated properties including tranmissivity, thermal
boundary conductance, and thermal conductivity
''' 
args = {'Gamma'}
      

def tau_spectral(Gamma, gb : AS, k, n_angle, T):
    '''
    Calculates a spectral tau which can be applied to the
    single mode, isotropic Callaway model. This relaxation 
    time is for kappa_xx.
    
    Don't think I need to adjust for the second grid. 
     
    '''
    d_angle = (math.pi / 2) / n_angle
    running_integrand = 0
    theta_list = np.arange(d_angle, math.pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, math.pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            running_integrand = running_integrand + \
            ((np.sin(theta - (d_angle / 2))**3 * np.cos(phi - (d_angle / 2))**2) * Gamma(k_vector_int)**(-1))
    
    return 8 * running_integrand * d_angle**2 * (3 / (4 * math.pi)) 
    

def transmissivity_spectral(Gamma, gb : AS, k, n_angle, T):        
    '''
    calculate the phonon spectral transmissivity
    '''
    vg = gb.vg_kmag(k)
    tau = tau_spectral(Gamma, gb, k, n_angle, T) * 1e-9
    a = vg * tau * gb.n_1D
    return a/((3/4)+a)


def Cv(k, T, omega_k):
    '''
    spectral heat capacity. there are no factors of group or phase velocity because it is in terms of k rather than omega
    But is it missing a factor of vg? maybe not if the integrand is dk?
    '''
    BEexp = hbar * omega_k/(kB*T)
    return (3*kB/(2*math.pi**2))*BEexp**2*k**2*math.exp(BEexp)/(math.exp(BEexp) - 1)**2 # This includes the k**2 from the volume integral.. very weird. I think I should change these.

   
def kL_spectral(Gamma, gb : AS, k, n_angle, T):
    '''
    calculate the spectral thermal conductivity
    '''
    vg = gb.vg_kmag(k)
    tau = tau_spectral(Gamma, gb, k, n_angle, T) * 1E-9
# "spectral heat capacity. There is no factors of group and phase velocity because it is in terms of k rather than omega.
    #print(Cv(k, T, omega_k))
    kL = (1/3) * Cv(k,T, gb.omega_kmag(k))*vg**2*tau
    return kL
    

def kL_T(Gamma, kmax, dk, vg_k, omega_k, T, n):
    '''
    Calculate the thermal conductivity versus temperature
    '''
    kLint = 0
    for k in np.arange(dk, kmax, dk):
        vg = vg_k(k)
        tau = tau_spectral(Gamma, k, vg_k, n) * 1E-9
        kLint = kLint + Cv(k,T, omega_k)*vg**2*tau*dk
    return (1/3)*kLint


def tbc_spectral(Gamma, gb : AS, k, n_angle, T):
    '''
    Calculate the spectral thermal boundary conductance

    '''
    vg = gb.vg_kmag(k)
    a = transmissivity_spectral(Gamma, gb, k, n_angle, T)
    tbc = a*vg*Cv(k, T, gb.omega_kmag(k))
    return (1/4)*tbc


def tbc_T(kmax, dk, vg_k, omega_k, T, n_1D, Gamma, n):
    '''
    Calculate the thermal boundary conductance versus T
    '''
    tbc_int = 0
    for k in np.arange(dk, kmax, dk):
        vg = vg_k(k)
        a = transmissivity_spectral(k, vg_k, n_1D, Gamma, n)
        tbc_int = tbc_int + a*vg*Cv(k,T, omega_k)*dk
    return (1/4)*tbc_int

def calculate_spectral_props(gb : AS, Gamma, prop_list = ['tau', 'transmissivity', 'TBC', 'kappa'],\
                             n_angle = 100, n_k = 100, T = 300):
    '''
    Calculate frequency-dependent properties
    prop_list : spectral properties which should be calculated
    '''
    spectral = {'vg' : [], 'omega' : []}
    function = {'tau' : tau_spectral, 'transmissivity' : transmissivity_spectral,
                'TBC' : tbc_spectral, 'kappa' : kL_spectral}
    if any(prop_list) not in ['tau', 'transmissivity', 'TBC', 'kappa']:
        ValueError('Property not in allowed values list')
    if 'tau' in prop_list:
        spectral['tau'] = []
    if 'transmissivity' in prop_list:
        spectral['transmissivity'] = []
    if 'TBC' in prop_list:
        spectral['TBC'] = []
    if 'kappa' in prop_list:
        spectral['kappa'] = []
    dk = gb.k_max / n_k
    k_mags = np.arange(dk, gb.k_max, dk)
    params = {'Gamma' : Gamma, 'gb' : gb, 'k' : dk, 'n_angle' : n_angle, 'T' : T}
    for k in k_mags: 
        spectral['vg'].append(gb.vg_kmag(k))
        spectral['omega'].append(gb.omega_kmag(k))
        params['k'] = k
        for prop in prop_list:
            spectral[prop].append(function[prop](**params))
    return spectral

