#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:09:04 2020

@author: ramyagurunathan

Helper functions
"""

from ThermalModel import Elasticity, CrystStruct, ThermalTransport
from math import pi as pi, cos, sin
import numpy as np

'''
Constants
'''
hbar = 6.6260696 * 10**-34 / (2 * pi)
kB = 1.38064852e-23

'''
Rotation Tensors
'''
rot_tensor_x = lambda theta: [[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]]
rot_tensor_y = lambda theta: [[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]]
rot_tensor_z = lambda theta: [[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0, 0, 1]]

def input_dict_from_MP(mpid):
    elast = Elasticity(mpid)
    struct = CrystStruct(mpid)
    thermal = ThermalTransport(struct, elast)
    input_dict = {'atmV' : struct.vol / struct.natoms,
                   'N' : struct.natoms,
                   'avg_vs' : thermal.therm_dict['average v_s'],
                   'gruneisen' : thermal.therm_dict['gruneisen'],
                   'nu' : elast.poisson,
                   }
    return input_dict

def k_mag(k_vector):
    '''
    Input: k_vector expressed as array or list of length 3
    Output: scalar magnitude value
    Description:
    Outputs magnitude of the k_vector
    '''
    return np.sqrt(k_vector[0]*k_vector[0] + k_vector[1]*k_vector[1] + k_vector[2]*k_vector[2])



def average_group_velocity(vg_mat):
   '''
   If passed a group velocity matrix of the acoustic phonons, provides the 
   average group velocity scalar value
   '''
   vt_3 = (vg_mat[0][0]**2 + vg_mat[0][1]**2 + vg_mat[0][2]**2)**(1/2)
   vt2_3 = (vg_mat[1][0]**2 + vg_mat[1][1]**2 + vg_mat[1][2]**2)**(1/2)
   vl_3 = (vg_mat[2][0]**2 + vg_mat[2][1]**2 + vg_mat[2][2]**2)**(1/2)
   avg_vg = ((1/3)*(1/vt_3**3 + 1/vt2_3**3 + 1/vl_3**3))**(-1/3)
   return avg_vg*1000

def average_phase_velocity(vp):
    avg_vp = ((1/3) * vp[0]**(-3) + vp[1]**(-3) + vp[2]**(-3))**(-1/3)
    return avg_vp*1000