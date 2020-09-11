#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 14:25:15 2020

@author: ramyagurunathan

Tilt Scattering Update AMM

"""

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import ScatteringPlots as SPlt
#from AngularVs import ElasticProps, v_long, v_shear, v_sound, v_xy
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
             'nu' : 0.29,
             'gruneisen' : 1,
        }

cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330

geom = 'tilt'
theta = 2

'''
Initialize input dictionary with Materials Project
'''
#input_dict = helper.input_dict_from_MP('mp-149')

#Instantiate as object of ArrayScattering
#tilt = AS(**input_dict, geom = 'tilt', theta = 5, ax = 1, d_GS = 350E-9)
#amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])

def initialize(input_dict, cmat, density, theta, geom, ax = 1, d_GS = 350e-9):
    amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])
    tilt = AS(**input_dict, geom = geom, amm = amm, theta = theta, ax = ax, d_GS = d_GS)
    return tilt


'''
Alternate rotation Scattering Potential
'''
#def V_R_ch_snell(k_vector, kprime_vector):
#    kmag = helper.k_mag(k_vector)
#    knorm = k_vector / kmag
#    v1, v2, theta2 = amm.vs_rot_Snell(knorm, helper.rot_tensor_z, tilt.theta)
#    if math.isnan(theta2):
#        qx = 2 * k_vector[0]
#    else:
#        qx = kmag * cos(theta2) - kmag * cos(acos(k_vector[0]/kmag)) 
#        if qx == 0:
#            return 0
#    print('rotate')
#    print(abs((helper.hbar / (math.sqrt(2 * pi))) * abs((v2 - v1)) * (kmag / qx))**2)
#    return  abs((helper.hbar / (math.sqrt(2 * pi))) * abs((v2 - v1)) * (kmag / qx))**2 #D is the width of the interface..??

'''
Rotation Scattering Potential
'''
def V_tilde_sq_R(tilt, k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
    v1, v2, theta2 = tilt.amm.vs_rot_Snell(knorm, helper.rot_tensor_z, tilt.theta)
#    if math.isnan(theta2):
#        qx = 2 * k_vector[0]
#    else:
#        qx = kmag * cos(theta2) - kmag * cos(acos(k_vector[0]/kmag))
#        if qx == 0:
#            return 0
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
# missing a negative sign from the factor of "i"?
    

def V1_twiddle_sq_S(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               (tilt.b / (1 - tilt.nu)) * ((q_vector[0] * q_vector[1]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2
               

def Gamma_GBS(tilt, k_vector, kprime_vectors):
#    print([tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax), tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)])
    return tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax) \
          + tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)\
          + tilt.GammaArray_rot(k_vector, V_tilde_sq_R)

def Gamma_GBS_rot_only(tilt, k_vector, kprime_yvectors, kprime_zvectors):
    return tilt.GammaArray_rot(k_vector, V_tilde_sq_R)

#Move to thermalTransport?
def Gamma(tilt, k_vector):
    return Gamma_GBS(tilt, k_vector, tilt.kprimes_y(k_vector)) * 1E-9 

def Gamma_rot_only(tilt, amm, k_vector):
    return Gamma_GBS_rot_only(tilt, k_vector, tilt.kprimes_y(k_vector), tilt.kprimes_z(k_vector)) * 1E-9

def calculate_Gammas(tilt, n_k):
    dk = tilt.k_max / n_k
    k_mags = np.arange(dk, tilt.k_max, dk)
    k_norm = k_mags / tilt.k_max
    k_vectors = []
    for k in k_mags:
        k_vectors.append([k, 0, 0])
    Gamma_GBS_list = []
    for k_vector in k_vectors:
        Gamma_GBS_list.append(Gamma(k_vector))
    return [k_norm, Gamma_GBS_list] 

if __name__ == "__main__":
    tilt = initialize(input_dict, cmat, density, theta, geom = 'tilt', ax = 1, d_GS = 350e-9)
#    Gamma_list = calculate_Gammas(tilt, 200)
#    SPlt.diffraction_plot(tilt, Gamma_list[0], Gamma_list[1])
    SPlt.convergence_tau_plot(tilt, Gamma, 110, T = 300, save = True)
#    spectral = TT.calculate_spectral_props(tilt, Gamma_rot_only, prop_list = ['tau'],\
#                                         n_angle = 200, n_k = 50, T = 300)
#    with open('spectral.json') as json_file:
#        spectral = json.load(json_file)
#    SPlt.spectral_plots(tilt, spectral, prop_list = ['tau'], save = True)

               
