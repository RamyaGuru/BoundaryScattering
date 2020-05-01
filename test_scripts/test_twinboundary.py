#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:57:53 2020

@author: ramyagurunathan

Twin Scattering Test file
"""

'''
Microstructural Features
'''
d_GS = 350E-9           # Average grain size [m]
n_1D = 3 / d_GS         # Number density of GBs [m^-1] is the 3 from dimensionality?
theta = 6.5 * (pi/180)

'''
Bi2Te3
'''
#cmat = np.array([[74.4, 21.7, 27.0, 13.3, 0.0, 0.0], [21.7, 74.4, 27.0, -13.3, 0.0, 0.0],\
#                [27.0, 27.0, 47.7, 0.0, 0.0, 0.0], [13.3, -13.3, 0.0, 27.4, 0.0, 0.0],\
#                [0.0, 0.0, 0.0, 0.0, 27.4, 13.3], [0.0, 0.0, 0.0, 0.0, 13.3, 26.4]])
#density = 7480

'''
Si
'''
cmat = np.array([[166.0, 64.0, 64.0, 0.0, 0.0, 0.0], [64.0, 166.0, 64.0, 0.0, 0.0, 0.0], [64.0, 64.0, 166.0, 0.0, 0.0, 0.0],\
                [0.0, 0.0, 0.0, 80.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 80.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 80.0]])

density = 2329

#Initialize a Christoffel object
ch = Christoffel(cmat, density)

#Set the z_direction to (0,0,1) and the x_direction to (1,0,0)
ch.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])

#Rotation tensordescribing the misorientation
rot_tensor_x = [[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]]
rot_tensor_y = [[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]]
rot_tensor_z = [[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0, 0, 1]]


'''
Input dict
'''
input_dict = {'avg_vs': average_group_velocity(vg_mat),
             'atmV': [2E-29],
             'N': 5,
             'nu' : 0.29,
             'gruneisen' : 1,
        }


