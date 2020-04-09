#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:32:44 2019

@author: ramyagurunathan

Plot strain fields of:
    1) tilt boundary array
    2) one of the twist boundary arrays

"""
from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

'''
Inputs
'''
vs = 6084.       # Speed of sound [m/s]
V = 2E-29       # Volume per atom [m^3]
N = 2           # Number of atoms per primitive unit cell

nu = 0.27       # Poisson's ratio

b = (V * N) ** (1 / 3)  # Burger's vector [m]

 


#Tilt Boundary real space strain fields

def epsDelta(x,y):
    return (-b / 2 * pi) * (1 - 2 * nu ) / (1 - nu) * (y / (x**2 + y**2))

def epsShear(x,y):
    return (b / 4 * pi) * x * (x**2 - y**2) / (x**2 + y**2)**2

def epsRot(x,y):
    return (b / pi) * (x / (x**2 + y**2))

    


#Twist boundary real space strain fields for the "m" - array

def eps13(x,y): 
    return -(b * y) / (4* pi * (x**2 + y**2))

def eps23(x, y):
    return (b * x) / (4 * pi * (x**2 + y**2))


#%% Main function for plotting the results for the twin boundary

nelem = 100
xlim = 1e-9
ylim = 1e-9
xrange = np.linspace(-xlim, xlim, 100)
yrange = np.linspace(-ylim, ylim, 100)

xval, yval = np.meshgrid(xrange, yrange)

strainDelta = epsDelta(xval, yval)
strainShear = epsShear(xval, yval)
strainRot = epsRot(xval, yval)

strain13 = eps13(xval, yval)
strain23 = eps23(xval, yval)
    
plt.pcolormesh(xval, yval, strain13, vmin = -.27, vmax = 0.27, cmap = 'Greys')
plt.colorbar()

#%%
'''
Apply appropriate periodicity to the atom displacements
'''

def u():
    
    


# Get the 