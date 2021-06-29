#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 14:25:15 2020

@author: ramyagurunathan

Symmetric Tilt Boundary Scattering 

"""

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import ScatteringPlots as SPlt
import math
import numpy as np
from math import asin, acos, pi, sin, cos, tan, atan
import helper
from AMMTransport import AMMTransport

np.seterr(divide='raise', invalid="raise")

'''
Input dict
'''
input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83,
             'shearmod' : 100,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330

geom = 'tilt'
theta = 5

'''
Initialize input dictionary with Materials Project
'''
#input_dict = helper.input_dict_from_MP('mp-149')

#Instantiate as object of ArrayScattering
#tilt = AS(**input_dict, geom = 'tilt', theta = 5, ax = 1, d_GS = 350E-9)
#amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])

'''
Initialize from input dictionary
'''

def initialize(input_dict, cmat, density, theta, geom, ax = 1, d_GS = 350e-9, bvK = True):
    amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])
    tilt = AS(**input_dict, geom = geom, amm = amm, theta = theta, ax = ax, d_GS = d_GS, bvK = bvK)
    return tilt


'''
Rotation Scattering Potential
'''
def V_tilde_sq_R(tilt, k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
    v1, v2, theta2 = tilt.amm.vs_rot_Snell(knorm, helper.rot_tensor_z, tilt.theta)
    return abs(helper.hbar * abs(v2 - v1) * abs(kmag / k_vector[0]))**2 * (2 * k_vector[0]**2 / kmag**2)

'''
Strain Field Scattering
'''

def V1_twiddle_sq_Delta(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               ((tilt.b * (1 - 2 * tilt.nu)) / (1 - tilt.nu)) * (q_vector[1]\
               / (q_vector[0]**2 + q_vector[1]**2)))**2
    

def V1_twiddle_sq_S(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               (tilt.b / (1 - tilt.nu)) * ((q_vector[0] * q_vector[1]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2
               
def V_core(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    return ((pi / 2) * 3 * helper.hbar * k *  tilt.V**(2/3))**2 * tilt.vs
               

def Gamma_GBS(tilt, k_vector, kprime_vectors):
#    print([tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax), tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)])
    return tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax) \
          + tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)\
          + tilt.GammaArray_rot(k_vector, V_tilde_sq_R)
          

def Gamma_GBS_rot_only(tilt, k_vector, kprime_yvectors, kprime_zvectors):
    return tilt.GammaArray_rot(k_vector, V_tilde_sq_R)

def Gamma_GBS_core(tilt, k_vector, kprime_vectors):
#    print([tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax), tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)])
    return tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax) \
          + tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)\
          + tilt.GammaArray(k_vector, kprime_vectors, V_core, tilt.ax)\
          + tilt.GammaArray_rot(k_vector, V_tilde_sq_R)
         
def Gamma_GBS_core_only(tilt, k_vector, kprime_vectors):
#    print([tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax), tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)])
    return tilt.GammaArray(k_vector, kprime_vectors, V_core, tilt.ax)\
          + tilt.GammaArray_rot(k_vector, V_tilde_sq_R)
          
          
def Gamma_GBS_core_no_rot(tilt, k_vector, kprime_vectors):
#    print([tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax), tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)])
    return tilt.GammaArray(k_vector, kprime_vectors, V_core, tilt.ax)

def Gamma(tilt, k_vector):
    return Gamma_GBS(tilt, k_vector, tilt.kprimes_y(k_vector)) * 1E-9 

def Gamma_rot_only(tilt, k_vector):
    return Gamma_GBS_rot_only(tilt, k_vector, tilt.kprimes_y(k_vector), tilt.kprimes_z(k_vector)) * 1E-9

def Gamma_core_only(tilt, k_vector):
    return Gamma_GBS_core_only(tilt, k_vector, tilt.kprimes_y(k_vector)) * 1E-9 

def Gamma_core(tilt, k_vector):
    return Gamma_GBS_core(tilt, k_vector, tilt.kprimes_y(k_vector)) * 1E-9 

def Gamma_core_no_rot(tilt, k_vector):
    return Gamma_GBS_core_no_rot(tilt, k_vector,  tilt.kprimes_y(k_vector)) * 1E-9



def calculate_Gammas(tilt, n_k, gamma_fxn = Gamma):
    dk = tilt.k_max / n_k
    k_mags = np.arange(dk, tilt.k_max, dk)
    k_norm = k_mags / tilt.k_max
    k_vectors = []
    for k in k_mags:
        k_vectors.append([k, 0, 0])
    Gamma_GBS_list = []
    for k_vector in k_vectors:
        Gamma_GBS_list.append(Gamma(tilt,k_vector))
    return [k_norm, Gamma_GBS_list] 

if __name__ == "__main__":
    theta = 12
    tilt = initialize(input_dict, cmat, density, theta, geom = 'tilt', ax = 1, d_GS = 350e-9)
    Gamma_list = calculate_Gammas(tilt, 200, Gamma_core_no_rot)
    SPlt.diffraction_plot(tilt, Gamma_list[0], Gamma_list[1])
#    SPlt.convergence_tau_plot(tilt, Gamma, 110, T = 300, save = True)
    spectral = TT.calculate_spectral_props(tilt, Gamma_core_no_rot, prop_list = ['tau'],\
                                         n_angle = 20, n_k = 10, T = 300, tau0 =  3.5)
#    with open('spectral.json') as json_file:
#        spectral = json.load(json_file)
    SPlt.spectral_plots(tilt, spectral, prop_list = ['tau'])

               
