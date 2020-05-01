#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:37:19 2019

@author: ramyagurunathan

Angular dependence of the phonon velocity on theta in a cubic material.

Required for estimation of acoustic mismatch due to anisotropy of the crystal
"""



'''
Elastic properties of the material
'''

from math import pi, cos, sin, tan
import numpy as np
import matplotlib.pyplot as plt
from ThermalModel import Elasticity, CrystStruct

        
'''
Functions: Speed of sound as a function of phi (shear and tranverse modes)
'''

'''
Initialize Elastic Properties from Materials Project
'''
#el = Elasticity('mp-149')
#struct = CrystStruct('mp-149')
#cmat = el.elastic_tensor
#density = struct.massDensity()

class ElastAnisotropy():
    def __init__(self, stiffness_matrix, density):
        self.cmat = stiffness_matrix
        self.rho = density
        self.C11 = self.cmat[0][0]
        self.C12 = self.cmat[0][1]
        self.C44 = self.cmat[3][3]

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
        
    '''
    Here, vx corresponds to the [010] direction and vy corresponds to the [001] direction
    '''
    
    def v_xy(self, pol, theta):
        vx = pol(theta)*cos(theta)
        vy = pol(theta)*sin(theta)
        return vx, vy

'''
Function: S parameter: change in speed of sound with a change in theta
'''

if __name__ == "__main__":
    pts = 200
    vt = np.zeros((pts, 2))
    vl= np.zeros((pts, 2))
    vs = np.zeros((pts, 2))
    '''
    Silicon Inputs
    '''
    stiffness_mat = [[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]]
    cmat = [[j * 1e9 for j in i] for i in stiffness_mat]
    density = 2330
    aniso = ElastAnisotropy(cmat, density)
    n=0
    for theta in np.linspace(0,2*pi, 200):
        vt[n][0], vt[n][1] = aniso.v_xy(aniso.v_shear_cubic, theta) 
        vl[n][0], vl[n][1] = aniso.v_xy(aniso.v_long_cubic, theta)
        vs[n][0], vs[n][1] = aniso.v_xy(aniso.v_sound_cubic, theta)
        n=n+1
    
    plt.plot(vl[:,0], vl[:,1], label = r'$v_l$', color = 'xkcd:blood red')
    plt.plot(vt[:,0], vt[:,1], label = r'$v_t$', color = 'xkcd:black')
    plt.plot(vs[:,0], vs[:,1], label = r'$v_{s, avg}$', color = 'xkcd:blue green')
    plt.xlim([-10000,10000])
    plt.ylim([-10000,10000])
    plt.ylabel(r'v$_y$ (m/s)')
    plt.xlabel(r'v$_x$ (m/s)')
    plt.legend(prop = {'size': 18})
    plt.savefig('angularvs.pdf', dpi=400, bbox_inches='tight')