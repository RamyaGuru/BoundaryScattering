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
import pickle


'''
Constants
'''
hbar = 6.6260696 * 10**-34 / (2 * math.pi)
kB = 1.38064852e-23


'''
PROPERTIES: options for calculated properties including tranmissivity, thermal
boundary conductance, and thermal conductivity
''' 
#      
def tau_spectral(Gamma, gb : AS, k, n_angle, T, directional = False):
    '''
    Calculates a spectral tau which can be applied to the
    single mode, isotropic Callaway model. This relaxation 
    time is for kappa_xx.
     
    '''
    d_angle = (math.pi / 2) / n_angle
    running_integrand = 0
    theta_list = np.arange(0, math.pi + d_angle, d_angle) # not including all incident angles? should check this
    phi_list = np.arange(0, math.pi + d_angle, d_angle)
    tau_directional = []
#    i = 0
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
#            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
#                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
#                            k * np.cos(theta - (d_angle / 2))]
            k_vector_int = [k * np.cos(theta - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2))]
            running_integrand = running_integrand + \
            ((2 * np.sin(theta - (d_angle / 2)) * np.cos(theta - (d_angle / 2))**2) * Gamma(k_vector_int)) * d_angle**2 # Integrate over scattering rate: See https://hackingmaterials.lbl.gov/amset/scattering/
            tau_directional.append([theta, phi, Gamma(k_vector_int)**(-1)])
#            i = i+1
#            if i == 10:
#                print(Gamma(k_vector_int))
#                i = 0
    if directional:
        with open('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/' + str(gb.geom) + str(gb.theta) + 'oldtau_directional_2_lf.pkl', 'wb') as file:
            pickle.dump(tau_directional, file)
    return ((3 / (4 * math.pi)) * running_integrand)**(-1) # need to think about this more.. the cos**2 now probably has to be inverted..
    

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
    kL = (1/3) * Cv(k,T, gb.omega_kmag(k))*vg**2*tau #because of 1/3.. can't I get rid of the kx/k term..? shouldn;t be there???
    return kL
    

def kL_T(Gamma, gb : AS, n_k, n_angle, T):
    '''
    Calculate the thermal conductivity versus temperature
    '''
    dk = gb.k_max/n_k
    kLint = 0
    for k in np.arange(dk, gb.k_max, dk):
        vg = gb.vg_kmag(k)
        tau = tau_spectral(Gamma, gb, k, n_angle, T) * 1E-9
        kLint = kLint + Cv(k,T, gb.omega_kmag(k))*vg**2*tau*dk
    return (1/3)*kLint


def tbc_spectral(Gamma, gb : AS, k, n_angle, T):
    '''
    Calculate the spectral thermal boundary conductance
    '''
    vg = gb.vg_kmag(k)
    a = transmissivity_spectral(Gamma, gb, k, n_angle, T)
    print(Cv(k, T, gb.omega_kmag(k)))
    tbc = (a / (1 - a))*vg*Cv(k, T, gb.omega_kmag(k))
    return (1/4)*tbc

def tbc_from_alpha(alpha, gb : AS, k, T):
    vg = gb.vg_kmag(k)
    return (1/4) *  (alpha / (1-alpha)) * vg * Cv(k, T, gb.omega_kmag(k))


def tbc_T(Gamma, gb : AS, n_k, n_angle, T):
    '''
    Calculate the thermal boundary conductance versus T
    '''
    dk = gb.k_max/n_k
    tbc_int = 0
    for k in np.arange(dk, gb.k_max, dk):
        vg = gb.vg_kmag(k)
        a = transmissivity_spectral(Gamma, gb, k, n_angle, T)
        tbc_int = tbc_int + (a / (1 - a))*vg*Cv(k,T, gb.omega_kmag(k))*dk
    return (1/4)*tbc_int



def calculate_spectral_props(gb : AS, Gamma, prop_list = ['tau', 'transmissivity', 'TBC', 'kappa'],\
                             function = {'tau' : tau_spectral, 'transmissivity' : transmissivity_spectral,
                'TBC' : tbc_spectral, 'kappa' : kL_spectral},
                             n_angle = 100, n_k = 100, T = 300):
    '''
    Calculate spectral properties
    prop_list : spectral properties which should be calculated
    
    Output:
        spectral: dictionary, key correspond to the property (vg and omega list are always included)
        value is the list of property values at the freqeuncies in the omega list
    '''
    spectral = {'vg' : [], 'omega' : []}
#    function = {'tau' : tau_spectral, 'transmissivity' : transmissivity_spectral,
#                'TBC' : tbc_spectral, 'kappa' : kL_spectral}
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

def calculate_temperature_dependence(gb : AS, Gamma, temp_list, prop_list = ['TBC', 'kappa'],\
                                     n_angle = 100, n_k = 50):
    '''
    Calculate properties at a range of temperatures
    
    Output:
        temp_dependence: dictionary, key corresponds to property (temperature always included as independent variable)
        value is the list of property values at different temperatures
    '''
    temp_dependence = {'temp' : temp_list}
    function = {'TBC': tbc_T, 'kappa': kL_T}
    if 'TBC' in prop_list:
        temp_dependence['TBC'] = []
    if 'kappa' in prop_list:
        temp_dependence['kappa'] = [] 
    params = {'Gamma' : Gamma, 'gb' : gb, 'n_k' : n_k, 'n_angle' : n_angle, 'T' : temp_list[0]}    
    for T in temp_list:
        params['T'] = T
        for prop in prop_list:
            temp_dependence[prop].append(function[prop](**params))
    return temp_dependence
    
def transport_coeffs_from_tau(gb : AS, k_list, tau_spectral, T, save = False):
    '''
    Return the thermal conductivity and TBC from the spectral tau array
    '''
    kappa = 0
    TBC = 0
    alpha = []
    dk =  gb.k_max / len(k_list)
    for tau, k in zip(tau_spectral, k_list):
        vg = gb.vg_kmag(k)
        Cv_s = Cv(k, T, gb.omega_kmag(k))
        kappa = kappa + Cv_s * vg**2 * tau * 1E-9 * dk
        a = (vg * tau * 1E-9 * gb.n_1D) / ((3/4) + (vg * tau * 1E-9 * gb.n_1D))
        alpha.append(a)
        TBC = TBC + (a / (1 - a)) * vg * Cv_s * dk
    transport = {'kappa': kappa / 3, 'TBC' : TBC / 4, 'spectral_alpha': alpha}
    if save:
        np.savez('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/datafiles/' +\
                 str(gb.geom) + str(gb.theta) + 'transport' + str(T) + '.npz')
    return transport

# add method to get the thermal boundary conductance from a scalar transmissivity value    

#def tbc_from_alpha(amm : AMM, alpha, n_k, T):
#    TBC = 0