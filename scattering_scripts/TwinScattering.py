#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 16:51:24 2020

@author: ramyagurunathan

Twin scattering script

Predct the transmissivity at a twin using the AMM Transmissivity treatment
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
             'nu' : 0.29,
             'gruneisen' : 1,
        }
'''
Silicon Data
'''
#cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
#     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])
#
#density = 2330

'''
Bismuth Relluride
'''
cmat = np.array([[74.4, 21.7, 27.0, 13.3, 0.0, 0.0], [21.7, 74.4, 27.0, -13.3, 0.0, 0.0],\
        [27.0, 27.0, 47.7, 0.0, 0.0, 0.0], [13.3, -13.3, 0.0, 27.4, 0.0, 0.0],\
        [0.0, 0.0, 0.0, 0.0, 27.4, 13.3], [0.0, 0.0, 0.0, 0.0, 13.3, 26.4]])

density = 7700

twin = AS(**input_dict, geom = 'twin', theta = 60, d_GS = 350E-9)
amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])

def V_R_ch_snell(k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
    v1, v2, theta2 = amm.vs_rot_Snell(knorm, helper.rot_tensor_z, twin.theta)
    if math.isnan(theta2):
        qx = 2 * k_vector[0]
    else:
        qx = kmag * cos(theta2)
    return  ((helper.hbar / (2 * pi)) * abs(1000 * (v2 - v1)) * (kmag / qx))**2      

def Gamma(k_vector):
    return twin.GammaArray_rot(k_vector, V_R_ch_snell) * 1E-9

def calculate_Gammas(n_k):
    dk = twin.k_max / n_k
    k_mags = np.arange(dk, twin.k_max, dk)
    k_norm = k_mags / twin.k_max
    k_vectors = []
    for k in k_mags:
        k_vectors.append([k, 0, 0])
    Gamma_GBS_list = []
    for k_vector in k_vectors:
        Gamma_GBS_list.append(Gamma(k_vector))
    return [k_norm, Gamma_GBS_list] 

if __name__ == "__main__":
    Gamma_list = calculate_Gammas(100)
    spectral = TT.calculate_spectral_props(twin, Gamma, prop_list = ['transmissivity', 'TBC'],\
                                         n_angle = 200, n_k = 20, T = 300)
    
