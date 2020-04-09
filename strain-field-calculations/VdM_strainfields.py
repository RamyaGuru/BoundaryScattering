#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 10:14:37 2018

@author: ramyaguru

Heterointerface stress fields: Van der Merwe model
"""
from math import pi as pi
from math import exp, sin, cos, tan, log
import matplotlib.pyplot as plt
import numpy as np

#%%
"""
Input Values for Silicon: From Riley's script-- maybe read from file later?
"""
a = 5.4e-10 #lattice parameter for Si 
a2 = 5.658e-10 #lattice parameter for Ge
b = 3.83e-10
nu = 0.27 #Poisson ratio
mu = 6.81e10 #shear modulus in Pascals
af = 1.56 # anisotropy factor 
kappa = 3/b
P = a/(a2-a) #characterizes misfit between lattices 1 and 2
p = P*b #spacing between misfit dislocations
c = ((P+1)*a)/(P+0.5) #displacement between lattices 1 and 2 over a distance p
#c= a
#p = 2e-9
beta = (pi*c)/(p*(1-nu)) # letting constant mu0=mu in Equation 13
#%%

"""
Function for Fourier coefficients
"""
def A(n):
    return ((1+ beta**2)**(1/2) - beta)**n
    
#%%
"""
Functions for the shear stress components-- from Van der Merwe
"""

def X(x):
    return (2*pi*x)/p

def Z(z):
    return (2*pi*z)/p

def R2(x, z, sgn):
    return 1+ A(2)*np.exp(sgn*2*Z(z)) - 2*A(1)*np.exp(sgn*Z(z))*np.cos(X(x))

def Q(x,z,sgn):
    return mu*c*A(1)*np.exp(sgn*Z(z))/((1-nu)*p*R2(x,z,sgn)**2)
    
def xhi(x, z, sgn):
    return ((mu*c*p)/(4*pi**2*(1-nu)))*Z(z)*np.log(R2(x,z,sgn)**(1/2))

def sigxx(x,z,sgn):
    return -Q(x,z,sgn)*(Z(z)*((1+A(2)*np.exp(sgn*2*Z(z)))*np.cos(X(x)) - 2*A(1)*np.exp(sgn*Z(z)))\
              +sgn*2*R2(x,z,sgn)*(np.cos(X(x)) - A(1)*np.exp(sgn*Z(z))))
def sigxx_riley(x,z,sgn):
    return Q(x,z,sgn)*Z(z)*((1+ A(2)*np.exp(sgn*2*Z(z)))*np.cos(X(x)) - 2*A(1)*np.exp(sgn*Z(z)))

def sigzz(x,z,sgn):
    return Q(x,z,sgn)*Z(z)*((1+ A(2)*np.exp(sgn*2*Z(z)))*np.cos(X(x)) - 2*A(1)*np.exp(sgn*Z(z)))

def sigzx(x,z,sgn):
    return -Q(x,z,sgn)*np.sin(X(x))*(R2(x,z,sgn) +sgn*Z(z)*(1-A(1)**2*np.exp(sgn*2*Z(z))))
#%%
"""
Function for the plane strain (exx)-- would combine dil, shear, and rot?
"""
#Defintion for the sigxx strain based on the Van der Merwe paper
def exx(x,z,sgn):
    return ((1-nu)*sigxx(x,z,sgn) - nu*sigzz(x,z,sgn))*(1/(2*mu))

#Definition for sigzz strain.. just flipped the dependences on the Poisson ratio 
def ezz(x,z,sgn):
    return ((1-nu)*sigzz(x,z,sgn) - nu*sigxx(x,z,sgn))*(1/(2*mu))

#Defintion of the shear strain.. just directly related to the shear stress for a cubic system
def ezx(x,z,sgn):
    return (1/mu)*sigzx(x,z,sgn)
#%%
"""
Function for the displacements
"""
    
#%%
"""
Main function: plotting the stress and strain fields
"""
xlim = 1
zlim = 1
xstrain = 1
zstrain = 1
nelem = 100
x = np.linspace(xlim* -1e-8,xlim*1e-8,nelem)
z = np.linspace(zlim *-1e-8,zlim* 1e-8,nelem)
Xax, Zax = np.meshgrid(x,z)

#Tighter window for strain
x_e= np.linspace(xstrain* -1e-8,xstrain*1e-8,nelem)
z_e= np.linspace(zstrain *-1e-8,zstrain* 1e-8,nelem)
Xstrain, Zstrain = np.meshgrid(x_e,z_e)


for s in [-1,1]:
    sigmaxx = np.zeros([0,nelem])
    sigmazx = np.zeros([0, nelem])
    sigmazz = np.zeros([0, nelem])
    strainxx = np.zeros([0, nelem])
    strainzz = np.zeros([0, nelem])
    strainzx = np.zeros([0, nelem])
    for zet in z: #Note: z is the row label and x is the column label
        sigmaxx_zet = sigxx(x,zet,s)
        sigmazx_zet = sigzx(x, zet,s)
        sigmazz_zet = sigzz(x,zet,s)
        strainxx_zet = exx(x, zet, s)
        strainzz_zet = ezz(x,zet,s)
        strainzx_zet = ezx(x,zet,s)
        sigmaxx_zet = np.transpose(sigmaxx_zet.reshape(nelem,1))
        sigmazx_zet = np.transpose(sigmaxx_zet.reshape(nelem,1))
        sigmazz_zet = np.transpose(sigmaxx_zet.reshape(nelem,1))
        strainxx_zet = np.transpose(strainxx_zet.reshape(nelem, 1))
        strainzz_zet = np.transpose(strainzz_zet.reshape(nelem, 1))
        strainzx_zet = np.transpose(strainzx_zet.reshape(nelem, 1))
        sigmaxx = np.append(sigmaxx,sigmaxx_zet, axis = 0)
        sigmazx = np.append(sigmazx, sigmazx_zet, axis = 0)
        sigmazz = np.append(sigmazz, sigmazz_zet, axis = 0)
        strainxx = np.append(strainxx, strainxx_zet, axis = 0)
        strainzz = np.append(strainzz, strainzz_zet, axis = 0)
        strainzx = np.append(strainzx, strainzx_zet, axis = 0)
    plt.figure()
    #v = np.linspace(-1e9, 1e9, 1000, endpoint=True)
    plt.pcolormesh(Xax, Zax, sigmaxx, vmin = -1e8, vmax = 1e8, cmap = "Greys")
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title(r'$\sigma_{xx}$')
    plt.figure()
    plt.pcolormesh(Xax, Zax, sigmazx, vmin = -1e8, vmax = 1e8, cmap = "Greys")
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title(r'$\sigma_{zx}$')
    plt.figure()
    plt.pcolormesh(Xax, Zax, sigmazz, vmin = -1e8, vmax = 1e8, cmap = "Greys")
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title(r'$\sigma_{zz}$')
    plt.figure()
    plt.pcolormesh(Xstrain, Zstrain, strainxx, vmin = 0, vmax = 0.27, cmap = "Greys")
    plt.ylim([-1e-9, 1e-9])
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title(r'$e_{xx}$')
    plt.colorbar()

#%% 2D plot of stress (sigma_xx) along the z-axis, fixed at x=0
plt.figure()    
plt.plot(z[:], sigmaxx[:,50])
plt.xlabel('z')
plt.ylabel(r'$\sigma_{xx}$')
plt.figure()
plt.plot(z[:50], sigmaxx[:50, 50])
plt.xlabel('z')
plt.ylabel(r'$\sigma_{xx}$')
plt.figure()
plt.xlabel('z')
plt.ylabel(r'$\sigma_{xx}$')
plt.plot(z[50:], sigmaxx[50:,50])

#%% 2D plot of strain (e_xx) along the z-axis, fixed at x=0
plt.figure()
plt.plot(z[25:75], strainxx[25:75, 50])
plt.xlabel('z')
plt.ylabel(r'$e_{xx}$')

plt.figure()
plt.plot(z[25:75], strainzz[25:75, 50])
plt.xlabel('z')
plt.ylabel(r'$e_{zz}$')

plt.figure()
plt.plot(z[25:75], strainzx[25:75, 50])
plt.xlabel('z')
plt.ylabel(r'$e_{zx}$')

#%% 2D plot of