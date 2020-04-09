#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:12:30 2019

@author: ramyagurunathan

Strain field plots for screw dislocation array of twist boundary
"""

from math import pi as pi
from math import exp, sin, cos, tan, log
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols
import sympy.plotting as smp

V = 2E-29 #Volume per atom
N = 2 #number of atoms in the unit cell
b = (V*N)**(1/3) #Burgers vector

#Define symbols for directions and n

x1, x2, x3, n, m = symbols('x1 x2 x3 n m')
#Functions to plot the strain fields

#Periodicity constraints
def X(x):
    return 2*pi*x/D

def Z(z):
    return 2*pi*z/D

#m-array strain fields, where the screw dislocations are pointed in the y-direction
def e12:
    return b*x3/(4*pi*(x1**2 + x3**2))
    
    
def e23:
    return b*x1/(4*pi*(x1**2 + x3**2))
    
    

#n-array strain fields, where the screw dislocations are pointed in the z-direction'
    
def e13:
    return -b*x2/(4*pi*(x1**2 + x2**2))
    
def e23:
    return -b*x1/(4*pi*(x1**2 + x2**2))
    
    
#Do the Fourier sums and plot
    

    