#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 13:14:30 2018

@author: guru

Script for implementing the Van der Merwe model for boundary scattering
at a heterointerface
"""

from sympy import symbols 
from sympy.plotting import plot
import sympy.plotting as smp
from math import pi as pi
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.titlesize']= 12

#%% Define the features of the dislocation array

#b = (4e-29)**(1/3) # Burger's vector in meters
b = 3.83e-10 #Burgers vector used in the VdM model 
nu = 0.27 # Poisson ratio
D = 8e-9 #value of p in the VdM model
hbar = 1.504e-34
#%%

#Define the symbolic variables for the 
x, y, n = symbols('x y n')

#Expression for dilatational strain  of a dislocation centered at y= nD
eps_n = (-b/ 2*pi)*((1-2*nu)/(1-nu))*((y-n*D)/(x**2 + (y-n*D)**2))

#Plot of dilatational strain field versus y
py = plot(eps_n.subs([(x,0), (n,0)]), (y, -1000,1000), ylim=[-1e-11, 1e-11], title = 'Dilatational strain for single dislocation')


#Plot of dilatational strain field versus x
px = plot(eps_n.subs([(y,1), (n,0)]), (x, -1000,1000), title = 'Dilatational strain for single dislocation')


#Plot of the dilatational strain versus x for an infinite sum over n
sumeps = 0*x
for m in range(-100,100):
    sumeps = sumeps + eps_n.subs([(x,1), (n,m)])
    
#Plot of the dilatational strain versus y for an infinite sum over n 
    
psum = plot(sumeps, (y,-100,100), title = 'Dilatational strain: finite sum over n')


#Combined strain field for all of the dislocations in the array. Strictly, this
# is a discrete sum over the dislocations in the array, but approximated as a real
# space integral over y from -inf to inf to yield function:

eps_tot =(-b/(2*D))*((1-2*nu)/(1-nu))*(x/abs(x))

pm1 = plot(eps_tot, (x,-20,20), title = 'Strain: Approx. sum over n as an integral')




