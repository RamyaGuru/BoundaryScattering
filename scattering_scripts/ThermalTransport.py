#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:37:01 2020

@author: ramyagurunathan

Transport Properties Module
"""

import numpy as np
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

def tau_spectral(Gamma, k, vg_func, n):
    '''
    Calculates a spectral tau which can be applied to the
    single mode, isotropic Callaway model. This relaxation 
    time is for kappa_xx.
    
    Don't think I need to adjust for the second grid. 
     
    '''
    d_angle = (math.pi / 2) / n
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
            ((np.sin(theta - (d_angle / 2))**3 * np.cos(phi - (d_angle / 2))**2) * Gamma(k_vector_int, vg_func(k_vector_int))**(-1))
    
    return 8 * running_integrand * d_angle**2 * (3 / (4 * math.pi)) 
    

def transmissivity(k, vg_k, n_1D, Gamma, n):        
    '''
    calculate the phonon spectral transmissivity
    '''
    vg = vg_k(k)
    tau = tau_spectral(Gamma, k, vg_k, n) * 1e-9
    a = vg*tau*n_1D
    return a/((3/4)+a)


def Cv(k, T, omega_k):
    '''
    spectral heat capacity. there are no factors of group or phase velocity because it is in terms of k rather than omega
    But is it missing a factor of vg? maybe not if the integrand is dk?
    '''
    k_vec = [k,0,0]
    BEexp = hbar*omega_k(k_vec)/(kB*T)
    return (3*kB/(2*pi**2))*BEexp**2*k**2*math.exp(BEexp)/(math.exp(BEexp) - 1)**2 # This include the k**2 from the volume integral.. very weird. I think I should change these.

   
def kL_spectral(Gamma, k, vg_k, omega_k, T, n):
    '''
    calculate the spectral thermal conductivity
    '''
    vg = vg_k(k)
    tau = tau_spectral(Gamma, k, vg_k, n) * 1E-9
# "spectral heat capacity. There is no factors of group and phase velocity because it is in terms of k rather than omega.
    #print(Cv(k, T, omega_k))
    kL = (1/3) * Cv(k,T, omega_k)*vg**2*tau
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


def tbc_spectral(k, vg_k, omega_k, T, Gamma, n_1D, n):
    '''
    Calculate the spectral thermal boundary conductance

    '''
    vg = vg_k(k)
    a = transmissivity(k, vg_k, n_1D, Gamma, n)
    tbc = a*vg*Cv(k, T, omega_k)
    return (1/4)*tbc


def tbc_T(kmax, dk, vg_k, omega_k, T, n_1D, Gamma, n):
    '''
    Calculate the thermal boundary conductance versus T
    '''
    tbc_int = 0
    for k in np.arange(dk, kmax, dk):
        vg = vg_k(k)
        a = transmissivity(k, vg_k, n_1D, Gamma, n)
        tbc_int = tbc_int + a*vg*Cv(k,T, omega_k)*dk
    return (1/4)*tbc_int

