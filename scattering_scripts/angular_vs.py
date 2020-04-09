#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:37:19 2019

@author: ramyagurunathan

Angular dependence of the phonon velocity on theta in a cubic material.
"""



'''
Elastic properties of the material
'''

from math import pi, cos, sin, tan
import numpy as np
import matplotlib.pyplot as plt

class ElasticProps:
    def __init__(self, props):
        [self.rho, self.C11, self.C12, self.C44] = props
        
'''
Functions: Speed of sound as a function of phi (shear and tranverse modes)
'''

p = ElasticProps([2329, 16.564e10, 6.394e10, 7.951e10])


def v_shear(theta):
    v_s = (2*p.rho)**(-1/2)*(p.C11 + p.C44 - ((p.C11 - p.C44)**2*cos(2*theta)**2 + (p.C12+p.C44)**2*sin(2*theta)**2)**(1/2))**(1/2)
    return v_s

def v_long(theta):
    v_l = (2*p.rho)**(-1/2)*(p.C11+p.C44 + ((p.C11-p.C44)**2*cos(2*theta)**2 + (p.C12 + p.C44)**2*sin(2*theta)**2)**(1/2))**(1/2)
    return v_l


'''
Here, vx corresponds to the [010] direction and vy corresponds to the [001] direction
'''

def v_xy(pol, theta):
    vx = pol(theta)*cos(theta)
    vy = pol(theta)*sin(theta)
    return vx, vy


pts = 200
vs = np.zeros((pts, 2))
vl= np.zeros((pts, 2))
n=0
for theta in np.linspace(0,2*pi, 200):
    vs[n][0], vs[n][1] = v_xy(v_shear, theta) 
    vl[n][0], vl[n][1] = v_xy(v_long, theta)
    n=n+1

plt.plot(vl[:,0], vl[:,1])
plt.plot(vs[:,0], vs[:,1])
plt.xlim([-10000,10000])
plt.ylim([-10000,10000])

