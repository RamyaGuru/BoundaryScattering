#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 14:30:44 2019

@author: ramyagurunathan

Het interface Scattering using VdM strain fields

Scattering due to a semicoherent interface
Should see diffraction peaks and should see the momentum conservation condition


Rather than using hte Van der Merwe equation, try just using the stress fields
of an edge dislocation?


"""


from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import ScatteringPlots as SPlt
import math
import numpy as np
from math import asin, acos, pi, sin, cos, tan, atan
import helper
from AMMTransport import HetAMMTransport as HA
import time

np.seterr(divide='raise', invalid="raise")

#Crystal Properties
'''
Input values for Silicon and Germanium
'''

input_dict = {
        'avg_vs': (6084 + 5400) / 2, #currently using average velocity between material 1 and material 2
             'atmV': [1.97E-29, 2.27E-29],
             'N': 2,
             'nu' : 0.29,
             'gruneisen' : 1,
        }


'''
Stiffness matrices for Si and Ge
'''
#Si
cmat1 = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])
density1 = 2330

#Ge
cmat2 = np.array([[126.0, 44.0, 44.0, 0, 0, 0],[44.0, 126.0, 44.0, 0, 0, 0],[44.0, 44.0, 126.0, 0, 0 ,0],\
     [0, 0, 0, 67.7, 0, 0],[0, 0, 0, 0, 67.7, 0],[0 ,0, 0, 0, 0, 67.7]])

density2 = 5323

geom = 'heterointerface'

'''
Iniitalize ArrayScattering and AMMTransport Objects
'''

#het = AS(**input_dict, geom = 'heterointerface', ax = {'n' : 1, 'm' : 2}, d_GS = 350E-9)
#amm = HA(cmat1, cmat2, density1, density2, christoffel = True)

def initialize(input_dict, cmat : list, density : list, geom, ax = {'n' : 1, 'm' : 2}, d_GS = 350e-9):
    amm = HA(cmat[0], cmat[1], density[0], density[1], christoffel = True)
    het = AS(**input_dict, geom = 'heterointerface', amm = amm, ax = ax, d_GS = d_GS)
    return het
    


'''
AMM Term
'''

def V_tilde_amm(het, k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
    v1, v2 = het.amm.delta_vs(knorm)
    return abs(helper.hbar * abs(v2 - v1) * abs(kmag / k_vector[0]))**2 * (2 * k_vector[0]**2 / kmag**2)


'''
Strain Scattering Functions
'''

# n array

def V_twiddle_sq_Delta_n(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               ((het.b * (1 - 2 * het.nu)) / (1 - het.nu)) * (q_vector[0]\
               / (q_vector[0]**2 + q_vector[1]**2)))**2
    

def V_twiddle_sq_S_n(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               (het.b / (1 - het.nu)) * ((q_vector[1] * q_vector[0]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2


def V_twiddle_sq_R_n(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               het.b * ((2 * q_vector[1]) / (q_vector[0]**2 + q_vector[1]**2)))**2
               
def V_twiddle_n(het, k_vector, kprime_vector):
    return V_twiddle_sq_Delta_n(het, k_vector, kprime_vector) + V_twiddle_sq_S_n(het, k_vector, kprime_vector) +\
V_twiddle_sq_R_n(het, k_vector, kprime_vector)


# m array

def V_twiddle_sq_Delta_m(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               ((het.b * (1 - 2 * het.nu)) / (1 - het.nu)) * (q_vector[0]\
               / (q_vector[0]**2 + q_vector[2]**2)))**2
    

def V_twiddle_sq_S_m(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               (het.b / (1 - het.nu)) * ((q_vector[2] * q_vector[0]**2)\
               / (q_vector[0]**2 + q_vector[2]**2)**2)) ** 2


def V_twiddle_sq_R_m(het, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * het.omega_kmag(k) * het.gruneisen * \
               het.b * ((2 * q_vector[2]) / (q_vector[0]**2 + q_vector[2]**2)))**2
               

def V_twiddle_m(het, k_vector, kprime_vector):
    return V_twiddle_sq_Delta_m(het, k_vector, kprime_vector) + V_twiddle_sq_S_m(het, k_vector, kprime_vector) +\
V_twiddle_sq_R_m(het, k_vector, kprime_vector)


#Integrate over the scattering potentials
def Gamma_GBS(het, k_vector, kprime_yvectors, kprime_zvectors):
   tot = [het.GammaArray(k_vector, kprime_yvectors, V_twiddle_n, het.ax['n']) \
          ,het.GammaArray(k_vector, kprime_zvectors, V_twiddle_m, het.ax['m']) \
          ,het.GammaArray_rot(k_vector, V_tilde_amm)]
   return sum(tot)

def Gamma_GBS_rot_only(het, k_vector, kprime_yvectors, kprime_zvectors):
    return het.GammaArray_rot(k_vector, V_tilde_amm)

def Gamma_rot(het, k_vector):
    return Gamma_GBS(het, k_vector, het.kprimes_y(k_vector), het.kprimes_z(k_vector)) * 1E-9

def Gamma_rot_only(het, k_vector):
    return Gamma_GBS_rot_only(het, k_vector, het.kprimes_y(k_vector), het.kprimes_z(k_vector)) * 1E-9

if __name__ == '__main__':
    het = initialize(input_dict, cmat = [cmat1, cmat2], density = [density1, density2],\
                     geom = geom)
    SPlt.convergence_tau_plot(het, Gamma_rot_only, 110, T = 300, save = True)
#    spectral = TT.calculate_spectral_props(het, Gamma_rot, prop_list = ['tau'],\
#                                         n_angle = 100, n_k = 100, T = 300)    



