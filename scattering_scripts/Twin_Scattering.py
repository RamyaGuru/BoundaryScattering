#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 09:50:47 2020

@author: ramyagurunathan

Coherent Twin Boundary

-> no dislocation array, just the misoritentation in elastic properties


This is a 60 degree twist boundary normal to the "x" direction
"""
import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/christoffel')

import ArrayScattering as AS
from christoffel import Christoffel 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from math import pi, cos, sin, acos, asin

np.seterr(divide='raise', invalid="raise")



'''
Crystal Properties--- Adjust for Bismuth Telluride
'''
hbar = 6.626e-34/(2*math.pi)
# Crystal properties
#vs = 6084.       # Speed of sound [m/s]
V = 3.40E-29       # Volume per atom [m^3]
N = 5         # Number of atoms per primitive unit cell
def gamma(k_vector):
    return 1           # Gruneissen parameter
nu = 0.27       # Poisson's ratio
k_max = (6 * pi**2 / (V * N))**(1 / 3)  # maximum k-vector (spherical BZ, only acoustic branches)
vs = 2000
omega_D = vs * k_max  # Debye frequency [Hz]


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

def average_group_velocity(vg_mat):
    vt_3 = (vg_mat[0][0]**2 + vg_mat[0][1]**2 + vg_mat[0][2]**2)**(3/2)
    vt2_3 = (vg_mat[1][0]**2 + vg_mat[1][1]**2 + vg_mat[1][2]**2)**(3/2)
    vl_3 = (vg_mat[2][0]**2 + vg_mat[2][1]**2 + vg_mat[2][2]**2)**(3/2)
    avg_vg = ((1/3)*(1/vt_3 + 1/vt2_3 + 1/vl_3))**(-1/3)
    return avg_vg*1000

def vs_k(k_vector):
    ch.set_direction_cartesian(k_vector)
    vs = ch.get_group_velocity()
    avg_vs = average_group_velocity(vs)
    return avg_vs

def omega_k(k_vector):
    if np.size(k_vector)>1:
        k = AS.k_mag(k_vector)
    else:
        k = k_vector
    return vs_k(k_vector) * k


def vs_rot(k_vector, theta):
    '''
    Method to calculate the average speed of sound change for a misorientation angle
    '''
    ch.set_direction_cartesian(k_vector)
    vs1 = ch.get_group_velocity()
    vs1 = average_group_velocity(vs1)
    ch.rotate_tensor(rot_tensor_y)
    ch.set_direction_cartesian(k_vector)
    vs2 = ch.get_group_velocity()
    vs2 = average_group_velocity(vs2)
    ch.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])
    return vs1, vs2

def AMM_transmissivity(k_vector, theta):
    '''
    Method to calculate the transmissivity from AMM. Then, can calculate the 
    '''
    k = AS.k_mag(k_vector)
    v1, v2 = vs_rot(k_vector, theta)
    theta1 = acos(k_vector[0]/k)
    try:
        theta2 = asin((v2/v1)*sin(theta1))
        a = (4*v1*v2*cos(theta1)*cos(theta2))/((v1*cos(theta1)\
                        + v2*cos(theta2))**2)
    except: 
        #must get reflection here.. so no transmissivity?
        a = 0
    return a

#def AMM_transmissivity(misorient, cos_inc):
#    v1, v2 = vs_rot(k_vector, misorient)
#    try:
#        cos_trans = asin((v2/v1) * sin(acos(cos_inc)))
#        a = (4 * v1 * v2 * cos_inc * cos_trans)/()
#    except:
#        a = 0
#    return a

def V_tilde_sq_R(k_vector, kprime_vector, theta):
    '''
    Method to calculate the scattering potential from the change in phonon velocity
    '''
#    k = AS.k_mag(k_vector)
#    q = np.asarray(kprime_vector)- np.asarray(k_vector)
    v1, v2 = vs_rot(k_vector, theta)
    return (hbar*abs(v1 - v2) * ( 1 / (2 * pi)))**2 #for now, assuming k/q_x is 1

'''
Tilt Boundary Scattering Rate
'''   
    
def tau_spectral(k, misorient, n):
    '''
    Need to take reciprocal of the Gamma
    '''
    d_angle = (pi / 2) / n
    running_integrand = 0
    theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            for theta2 in theta_list:
                for phi2 in phi_list:
                    k_vector_prime = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))] # should have same magnitude of k
                    running_integrand = running_integrand  +\
                    V_tilde_sq_R(k_vector_int, k_vector_prime, misorient) * sin(theta - (d_angle/2)) * sin(theta2 - (d_angle/2))
    return 8 * running_integrand * d_angle**4 * 1e-9 * (3/ (4 * pi)) #Likely need to multiply by some coefficients? Riley did 8 * (3/4pi)

#                    (sin(theta2 - (d_angle/2)) * cos(phi2 - (d_angle/2)) - sin(theta1 - (d_angle/2)) * cos(phi2 - (d_angle/2)))**(-2)


def spectral_thermal_conductance(k, vs_k, omega_k, misorient, T, n):
    d_angle = (pi / 2) / n
    running_integrand = 0
    theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            a = AMM_transmissivity(k_vector_int, misorient)
            running_integrand = running_integrand + \
            sin(theta - (d_angle/2))**2 * cos(phi - (d_angle / 2)) * a * vs_k(k_vector_int)
    print(AS.Cv(k,T,omega_k))
    return (1/4) * AS.Cv(k, T, omega_k) * running_integrand * d_angle**2

def spectral_transmissivity(k, vs_k, omega_k, misorient, T, n):
    d_angle = (pi / 2) / n
    running_integrand = 0
    theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            a = AMM_transmissivity(k_vector_int, misorient)
            running_integrand = running_integrand + \
            (sin(theta - (d_angle / 2))**(2) * cos(phi - (d_angle / 2))) * a #Double check these powers? Shouldn't actually be squared?
    return running_integrand * d_angle**2

def spectral_kL(k, vs_k, omega_k, misorient, T, n):
    d_angle = (pi / 2) / n
    running_integrand = 0
    theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            a = AMM_transmissivity(k_vector_int, misorient)
            running_integrand = running_integrand + \
            sin(theta - (d_angle/2))**2 * cos(phi - (d_angle / 2)) * a * vs_k(k_vector_int)**2
    print(AS.Cv(k,T,omega_k))
    return (1/3) * AS.Cv(k, T, omega_k) * running_integrand * d_angle**2

'''
Initialize dictionary of spectral values
'''

    
def spectral_values(k, vs_k, omega_k, misorient, T, n_1D, n):
    '''
    Computes all spectral values startign from the transmissivity
    '''
    spectral_props = {
        'transmissivity' : [],
        'thermal boundary conductance' : [],
        'relaxation time' : [],
        'thermal conductivity': []}
    d_angle = (pi / 2) / n
    running_integrand = 0
    tbc_integrand = 0
    tcond_integrand = 0
    tau_integrand = 0
    theta_list = np.arange(d_angle, pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            a = AMM_transmissivity(k_vector_int, misorient)
            running_integrand = running_integrand + \
            (sin(theta - (d_angle / 2))**(2) * cos(phi - (d_angle / 2))) * a  
            tbc_integrand = tbc_integrand +\
            (sin(theta - (d_angle / 2))**(2) * cos(phi - (d_angle / 2))) * a * vs_k(k_vector_int)
            tau_integrand = tau_integrand + (3/4) * a * \
            (1/(vs_k(k_vector_int) * n_1D * (1 - a))) * vs_k(k_vector_int)**2
    spectral_props['transmissivity'] = running_integrand * d_angle**2
    spectral_props['TBC'] = (1/4) * AS.Cv(k, T, omega_k) * tbc_integrand * d_angle**2 
    spectral_props['tau'] = (1/3) * AS.Cv(k, T, omega_k) * tau_integrand
    return spectral_props

'''
Main function
'''

n_k = 3
dk = k_max / n_k
k_mags = np.arange(dk, k_max, dk)
k_norm = k_mags / k_max
k_vectors = []

omega_list = []
vg_list = []
tbc_list = []
tau_list = []
kL_list = []
trans_list = []
TBC_list = []

T = 300
for k in k_mags:
    omega_list.append(omega_k([k,0,0])) # omega and vg are supposed to be only a function of k, not k_vector. This is tacky and needs to be fixed!
    vg_list.append(vs_k([k,0,0]))
#    tbc_list.append(spectral_thermal_conductance(k, vs_k, omega_k, theta, T, 50))
#    tau_list.append(tau_spectral(k, theta, 10)) 
    a = spectral_transmissivity(k, vs_k, omega_k, theta, T, 50)
    trans_list.append(a)
    TBC_list.append(spectral_thermal_conductance(k, vs_k, omega_k, theta, T, 50))
    
    
    
def plot_spectral_props(spectral_props: dict):
    for k,v in spectral_props:
        plt.figure()
        plt.scatter()
        plt.xlabel(r'\omega')
        plt.ylabel(k)
        

#%%  Thermal boundary Conductance Plot
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
#plt.ylabel(r'$TBC$')
#plt.plot(k_mags, tbc_list)
#plt.savefig('tiltBoundary_tbc.pdf', dpi=400, bbox_inches = 'tight')
#plt.show()            
#%% Relaxation Time Plot
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
#plt.ylabel(r'$\tau$; (s)')
#plt.plot(k_mags, tau_list)
#plt.savefig('tiltBoundary_tau.pdf', bbox_inches = 'tight')
#plt.show()

#%% Transmissivity Plot
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
#plt.ylabel('Transmissivity')
#plt.plot(k_mags, trans_list)
#plt.savefig('twin_trans.pdf', dpi=400, bbox_inches = 'tight')
#plt.show() 
#    
##%% Thermal boundary Conductance
#    
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
#plt.ylabel('Thermal Boundary Conductance')
#plt.plot(k_mags, TBC_list)
#plt.savefig('twin_TBC.pdf', dpi=400, bbox_inches = 'tight')
#plt.show() 