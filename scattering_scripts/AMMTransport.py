#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:00:55 2020

@author: ramyagurunathan

AMM Transport: Thermal Resistance due to Acoustic Mismatch at a boundary
because of the solid rotation and elastic anisotropy of the material
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/christoffel')

from christoffel import Christoffel 
import helper as h
import numpy as np
from math import cos, sin, acos, asin, pi, nan
import ThermalTransport as TT

class AMMTransport():
    def __init__(self, stiffness_matrix, density, atmV, N, christoffel = True):
        self.cmat = stiffness_matrix
        self.density = density
        self.C11 = self.cmat[0][0]
        self.C12 = self.cmat[0][1]
        self.C44 = self.cmat[3][3] 
#        self.theta = theta * (pi / 180)
        self.k_max = (6 * pi**2 / (atmV * N))**(1 / 3)
        if christoffel:
            self.ch_obj = Christoffel(stiffness_matrix, density)
            self.ch_obj.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])
        
    '''
    Simplified methods for Cubic Systems
    '''
    def v_shear_cubic(self, theta):
        v_t = (2*self.rho)**(-1/2)*(self.C11 + self.C44 - ((self.C11 - self.C44)**2*cos(2*theta)**2\
                             + (self.C12 + self.C44)**2*sin(2*theta)**2)**(1/2))**(1/2)
        return v_t
    
    def v_long_cubic(self, theta):
        v_l = (2*self.rho)**(-1/2)*(self.C11+self.C44 + ((self.C11-self.C44)**2*cos(2*theta)**2\
                          + (self.C12 + self.C44)**2*sin(2*theta)**2)**(1/2))**(1/2)
        return v_l
    
    #try the cubic mean thing instead from density of states
    def v_sound_cubic(self, theta):
        v_l = self.v_long(theta)
        v_t = self.v_shear(theta)
        v_s = ((1 / (3 * v_l**3)) + (2 / (3 * v_t**3)))**(-1/3)
        return v_s
        
    
    def v_xy(self, pol, theta):
        '''
        For generation of polar slowness plots. Here, vx corresponds to the [001] direction
        and vy corresponds to the [010] dorection
        '''
        vx = pol(theta)*cos(theta)
        vy = pol(theta)*sin(theta)
        return vx, vy
    
    '''
    General Methods for Arbitrary Symmetry
    '''
   
    def vs_k(self, k_vector):
        self.ch_obj.set_direction_cartesian(k_vector)
        vs = self.ch_obj.get_group_velocity()
        avg_vs = h.average_group_velocity(vs)
        return avg_vs
    
    def omega_k(self, k_vector):
        if np.size(k_vector)>1:
            k = h.k_mag(k_vector)
        else:
            k = k_vector
        return self.vs_k(k_vector) * k
    
    def vs_rot(self, k_vector, rot_tensor_method, theta):
        '''
        Method to calculate the average speed of sound change for a misorientation angle
        '''
        self.ch_obj.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])
        k_norm = k_vector / h.k_mag(k_vector)
        self.ch_obj.rotate_tensor(rot_tensor_method(-1 * (theta / 2)))
        self.ch_obj.set_direction_cartesian(k_norm)
        vs1 = self.ch_obj.get_group_velocity()
        vs1 = h.average_group_velocity(vs1) #a
        self.ch_obj.rotate_tensor(rot_tensor_method(theta / 2))
        self.ch_obj.set_direction_cartesian(k_norm)
        vs2 = self.ch_obj.get_group_velocity()
        vs2 = h.average_group_velocity(vs2)
        self.ch_obj.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])
        return vs1, vs2
    
    def vs_rot_fcc(self, k_vector, rot_tensor_method, theta):
        self.ch_obj.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [1.0, 1.0, 1.0])
        k_norm = k_vector / h.k_mag(k_vector)
        self.ch_obj.rotate_tensor(rot_tensor_method(-1 * (theta / 2)))
        self.ch_obj.set_direction_cartesian(k_norm)
        vs1 = self.ch_obj.get_group_velocity()
        vs1 = h.average_group_velocity(vs1) #a
        self.ch_obj.rotate_tensor(rot_tensor_method(theta / 2))
        self.ch_obj.set_direction_cartesian(k_norm)
        vs2 = self.ch_obj.get_group_velocity()
        vs2 = h.average_group_velocity(vs2)
        self.ch_obj.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [1.0, 1.0, 1.0])        
        return vs1, vs2

        
    
    def vs_rot_Snell(self, k_vector, rot_tensor_method, theta):
        k = h.k_mag(k_vector)
        v1, v2 = self.vs_rot(k_vector, h.rot_tensor_y, theta)
        theta1 = acos(k_vector[0]/k)
        try:
            theta2 = asin((v2/v1)*sin(theta1))
        except:
            theta2 = nan
        return v1, v2, theta2
    
    def AMM_transmissivity(self, k_vector, rotate_tensor_method, theta):
        '''
        Method to calculate the transmissivity from AMM. 
        '''
        k = h.k_mag(k_vector)
        v1, v2 = self.vs_rot(k_vector, rotate_tensor_method, theta)
        theta1 = acos(k_vector[0]/k)
        try:
            theta2 = asin((v2/v1)*sin(theta1))
            a = (4 * (v2 / v1) * (cos(theta2) / cos(theta1))) / ((v2/v1) + cos(theta2) / cos(theta1))**2
#            a = (4*v1*v2*cos(theta1)*cos(theta2))/((v1*cos(theta1)\
#                            + v2*cos(theta2))**2)
        except: 
            #must get reflection here.. so no transmissivity?
            a = 0
        return a
    
#    def AMM_transmissivity_fcc(self, ):
        
        
    
    def spectral_transmissivity(self, k, theta, rotate_tensor_method, T, n_angle = 100):
        d_angle = (pi / 2) / n_angle
        running_integrand = 0
        theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
        phi_list = np.arange(d_angle, 2 * pi, d_angle)
        for theta in theta_list:
            for phi in phi_list:
                #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
                k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                                k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                                k * np.cos(theta - (d_angle / 2))]
                a = self.AMM_transmissivity(k_vector_int, rotate_tensor_method, theta)
                running_integrand = running_integrand + \
                (sin(theta - (d_angle / 2))**(2) * cos(phi - (d_angle / 2))) * a #Double check these powers? Shouldn't actually be squared?.. a forward scattering suprresion thing?
        return running_integrand * d_angle**2
    
    def calculate_spectral_props(self, theta, rot_tensor_method, prop_list = ['transmissivity'],\
                                 function = {'transmissivity' : spectral_transmissivity},
                                 n_angle = 100, n_k = 100, T = 300):
        spectral = {'vg' : [], 'omega' : []}
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
        dk = self.k_max / n_k
        k_mags = np.arange(dk, self.k_max, dk)
        params = {'k' : dk, 'theta' : theta * (pi / 180), 'rotate_tensor_method' : rot_tensor_method, 'T' : T, 'n_angle' : n_angle, }
        for k in k_mags: 
            spectral['vg'].append(self.vs_k([k,0,0]))
            spectral['omega'].append(self.omega_k([k,0,0]))
            params['k'] = k
            for prop in prop_list:
                spectral[prop].append(function[prop](self, **params))
        return spectral
    
    def tbc_from_alpha(self, vs, alpha, n_k, T):
        TBC = 0
        dk = self.k_max / n_k
        k_mags = np.arange(dk, self.k_max, dk)
        for k in k_mags:
            TBC = TBC + TT.Cv(k, T, vs * k) * dk
        return (1/4) * (alpha/(1-alpha)) * vs * TBC
        
            
                
        
if __name__ == "__main__":
    Bi2Te3 = np.array([[74.4, 21.7, 27.0, 13.3, 0.0, 0.0], [21.7, 74.4, 27.0, -13.3, 0.0, 0.0],\
            [27.0, 27.0, 47.7, 0.0, 0.0, 0.0], [13.3, -13.3, 0.0, 27.4, 0.0, 0.0],\
            [0.0, 0.0, 0.0, 0.0, 27.4, 13.3], [0.0, 0.0, 0.0, 0.0, 13.3, 26.4]])
    Si = np.array([[166.0, 64.0, 64.0, 0.0, 0.0, 0.0], [64.0, 166.0, 64.0, 0.0, 0.0, 0.0], [64.0, 64.0, 166.0, 0.0, 0.0, 0.0],\
                [0.0, 0.0, 0.0, 80.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 80.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 80.0]])
    diamond = np.array([[1026.2, 90.558, 90.558, 0.0, 0.0, 0.0], [90.558, 1026.2, 90.558, 0.0, 0.0, 0.0], [90.558, 90.558, 1026.2, 0.0, 0.0, 0.0],\
                [0.0, 0.0, 0.0, 314.3, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 314.3, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 314.3]])
    InP = np.array([[87, 46, 46, 0.0, 0.0, 0.0], [46, 87, 46, 0.0, 0.0, 0.0], [46, 46, 87, 0.0, 0.0, 0.0],\
                [0.0, 0.0, 0.0, 42, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0,42, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 42]])
                    
                    
    '''
    Silicon
    '''
#    cmat = Si
#    density = 2329
#    N = 2
#    atmV = 40.89E-30 / 2
#    vs = 6084
    '''
    InP
    '''
#    cmat = InP
#    density = 4810
#    N = 2
#    atmV = 52.840E-30 / 2
#    vs = 3212
    '''
    Bi2Te3
    '''
    cmat = Bi2Te3
    density  = 7700
    N = 5
    atmV = 35.6E-30
    vs = 1767
    '''
    Diamond
    '''
#    cmat = diamond
#    density = 3510
#    N = 2
#    atmV = 11.41E-30 / 2
#    vs = 10542
    twin = AMMTransport(cmat, density, atmV = atmV, N = N, christoffel = True)
    spectral = twin.calculate_spectral_props(60, h.rot_tensor_z, prop_list = ['transmissivity'], n_angle = 300, n_k = 2, T = 300)
    TBC = twin.tbc_from_alpha(vs, spectral['transmissivity'][0], 500, 300)