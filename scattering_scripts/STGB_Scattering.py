#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:44:05 2019

@author: ramyagurunathan

STGB scattering

Steps:
    1) Set up a step function in the misorientation angle theta
    2) Set up the strain field 
    3) 
    
    Si-Si twist boundary
"""


import ArrayScattering as AS
from AngularVs import ElasticProps, v_long, v_shear, v_sound, v_xy
import Callaway as Cal
import math
import numpy as np
import matplotlib.pyplot as plt
from math import asin, acos, pi, sin, cos, tan, atan

np.seterr(divide='raise', invalid="raise")


'''
Crystal Properties
'''
hbar = 6.626e-34/(2*math.pi)
# Crystal properties
vs = 6084.       # Speed of sound [m/s]
V = 2E-29       # Volume per atom [m^3]
N = 2           # Number of atoms per primitive unit cell
def gamma(k_vector):
    return 1           # Gruneissen parameter
nu = 0.27       # Poisson's ratio
k_max = (6 * math.pi**2 / (V * N))**(1 / 3)  # maximum k-vector (spherical BZ, only acoustic brances)
omega_D = vs * k_max  # Debye frequency [Hz]
def omega_k(k_vector):
    if np.size(k_vector)>1:
        k = AS.k_mag(k_vector)
    else:
        k = k_vector
    return vs * k
def vg_k(k_vector):
    return vs

cmat = [[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]]

density = 2330

'''
Microstructural Features
'''

## Microstructural parameters
b = (V * N) ** (1 / 3)  # Burger's vector [m]
d_GS = 350E-9           # Average grain size [m]
n_1D = 3 / d_GS         # Number density of GBs [m^-1] is the 3 from dimensionality?
#D = 1E-9                # Linear defects spacing [m]

#Defect Spacing from the misorientation angle
twist = 70.4
theta =  twist * (pi / 180)#Misorientation angle

if theta:
    D = b / ( 2 * tan(theta / 2))
elif D:
    theta = 2*atan(b/(2*D))

'''
Angular dependence of the speed of sound 
'''
#ch = Christoffel(cmat, density)
#
##Set the z_direction to (0,0,1) and the x_direction to (1,0,0)
#ch.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])

def average_group_velocity(vg_mat):
    vt_3 = (vg_mat[0][0]**2 + vg_mat[0][1]**2 + vg_mat[0][2]**2)**(3/2)
    vt2_3 = (vg_mat[1][0]**2 + vg_mat[1][1]**2 + vg_mat[1][2]**2)**(3/2)
    vl_3 = (vg_mat[2][0]**2 + vg_mat[2][1]**2 + vg_mat[2][2]**2)**(3/2)
    avg_vg = ((1/3)*(1/vt_3 + 1/vt2_3 + 1/vl_3))**(-1/3)
    return avg_vg*1000

#def vs_k(k_vector):
#    ch.set_direction_cartesian(k_vector)
#    vs = ch.get_group_velocity()
#    avg_vs = average_group_velocity(vs)
#    return avg_vs*1000
#
#def vs_rot(k_vector, theta):
#    '''
#    Method to calculate the average speed of sound change for a misorientation angle
#    '''
#    ch.set_direction_cartesian(k_vector)
#    vs1 = ch.get_group_velocity()
#    vs1 = average_group_velocity(vs1)
#    ch.rotate_tensor([[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0, 0, 1]])
#    ch.set_direction_cartesian(k_vector)
#    vs2 = ch.get_group_velocity()
#    vs2 = average_group_velocity(vs2)
#    ch.rotate_tensor(x_dir = [1.0, 0.0, 0.0], z_dir = [0.0, 0.0, 1.0])
#    return vs1*1000, vs2*1000

'''
Rotation component:
'''

'''
theta here is the misorientation angle WHERE IS THETA??
inc: angle of incidence determined from the k-vector
direction normal to the plane is the x-axis


Change in velocity should come from misorientation angle rather than angle of incidence

Modify to do christoffel calculation
'''
def V_twiddle_sq_R(k_vector):
#    k = AS.k_mag(k_vector)
#    q_vector = np.asarray(kprime_vector)- np.asarray(k_vector)
#    #calculate the angle of incidence, what about misorientation angle?
#    inc = acos(k_vector[0]/(k))
    vs_new = v_sound(theta)
    # k will cancel out
    return (hbar*abs(vs_new - vs)/(2*math.pi))**2 #multiply by the ratio of k to q?

#def V_R_christoffel(k_vector, kprime_vector):
#    k = AS.k_mag(k_vector)
#    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
#    v1, v2 = vs_rot(k_vector, theta)
#    return (hbar / (2 * pi) * abs(v2 - v1) * (k / q_vector[0]))**2
    

'''
Strain components of the "n" array: Dislocaiton line in z, spacing in y
'''

def Vn_twiddle_sq_E13(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[1])/(2* (q_vector[0]**2 + q_vector[1]**2)))**2
    
def Vn_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[0]/(2*(q_vector[0]**2 + q_vector[1]**2))))**2

def V_twiddle_sq_n(k_vector, kprime_vector):
    return Vn_twiddle_sq_E13(k_vector, kprime_vector) + Vn_twiddle_sq_E23(k_vector, kprime_vector)
    

'''
Strain components of the "m" array: Dislocation line in y, spacing in z
'''
def Vm_twiddle_sq_E12(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[2])/(2* (q_vector[0]**2 + q_vector[2]**2)))**2
    
def Vm_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[0]/(2*(q_vector[0]**2 + q_vector[2]**2))))**2

def V_twiddle_sq_m(k_vector, kprime_vector):
    return Vm_twiddle_sq_E12(k_vector, kprime_vector) + Vm_twiddle_sq_E23(k_vector, kprime_vector)

'''
STGB scattering rates
'''
#Still need to include the rotation term
def Gamma_GBS(k_vector, kprime_yvectors, kprime_zvectors, vg, n_1D, D):
   return AS.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, vg, n_1D, D, 1) \
          + AS.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, vg, n_1D, D, 2) 


def Gamma_GBS_rot(k_vector, kprime_yvectors, kprime_zvectors, vg, n_1D, D):
   return AS.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, vg, n_1D, D, 1) \
          + AS.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, vg, n_1D, D, 2) \
          + V_twiddle_sq_R(k_vector)
          
def Gamma(k_vector, vg):
    return Gamma_GBS(k_vector, AS.kprimes_y(k_vector, D), AS.kprimes_z(k_vector, D), vg, n_1D, D) * 1E-9 #what's this function for?

def Gamma_rot(k_vector, vg):
    return Gamma_GBS_rot(k_vector, AS.kprimes_y(k_vector, D), AS.kprimes_z(k_vector, D), vg, n_1D, D) * 1E-9 
#NEED TO MODIFY
k_vector = [0.2 * k_max, 0, 0]

'''
DIFFRACTION CONDITION PLOT:
 Plot Gamma_GBS(k_vector) for normal incidence
'''
#%%
n_k = 1.E2
dk = k_max / n_k
k_mags = np.arange(dk, k_max, dk)
k_norm = k_mags / k_max
k_vectors = []
#for k in k_mags:
#    k_vectors.append([k, 0, 0]) #the non-zero term in this could be a problem? perpendicular to axis?
#
#Gamma_GBS_list = []
#for k_vector in k_vectors:
#    #Gamma_GBS_list.append(Gamma(k_vector, vg_k(k_vector)))
#    Gamma_GBS_list.append(Gamma_rot(k_vector, vg_k(k_vector)))
#        
#plt.figure()
#plt.xlim((0, 1))
#plt.ylim((0, 20))
#plt.xlabel(r'$k/k_{\mathrm{max}}$', fontsize=16)
#plt.ylabel(r'$\Gamma \; \mathrm{(ns^{-1})}$', fontsize=16)
#plt.plot(k_norm, Gamma_GBS_list)
#plt.savefig('twist_diffD1e-09.pdf', dpi=400, bbox_inches='tight')
#plt.show(block=False)
#%%
# Convergence of tau_spectral, n_angle=100 is sufficient.
# n_angle_list = np.arange(4, 100, 2)
# tau_nlist = []
# for n_angle in n_angle_list:
#     tau_nlist.append(AS.tau_spectral(Gamma, k_max / 5., vg_k, n_angle))

# plt.figure()
# plt.xlabel('n', fontsize=16)
# plt.ylabel(r'$\tau(k)^{-1} \; \mathrm{(ns^{-1})}$', fontsize=16)
# plt.plot(n_angle_list, tau_nlist)
# plt.show(block=False)

#%%
'''
Calculation of spectral tau and kappa


'''
#omega_list = []
#vg_list = []
#tau_list = []
#kappa_list = []
#trans_list = []
#tbc_list = []
#T = 300
#for k in k_mags:
#    omega_list.append(omega_k([k,0,0])) # omega and vg are supposed to be only a function of k, not k_vector. This is tacky and needs to be fixed!
#    vg_list.append(vg_k([k,0,0]))
#    tau_list.append(AS.tau_spectral(Gamma_rot, k, vg_k, 50))
#    trans = AS.transmissivity(k, vg_k, n_1D, Gamma_rot, 50)
#    trans_list.append(trans)
#    tbc_list.append(AS.tbc_spectral(k, vg_k, omega_k, T, Gamma_rot, n_1D, 50))
#    kappa_list.append(AS.kL_spectral(Gamma_rot, k, vg_k, omega_k, T, 50))
#    
#outfile = str(twist) + '_twist_vsf.npz'
#np.savez(outfile,np.array(tau_list), np.array(kappa_list), np.array(trans_list), np.array(tbc_list))
###
###
##    
###Save the generated lists
##outfile = str(twist) + '_spectral.npz'
##np.savez(np.array(omega_list), np.array(tau_list), np.array(tbc_list), np.array(kappa_list))
##    
##
##plt.figure()
##plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
##plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)
##plt.plot(k_mags, tau_list)
##plt.savefig('twistBoundary_D1e-9.pdf', dpi=400, bbox_inches = 'tight')
##plt.show(block=False)
##
##plt.figure
##plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
##plt.ylabel(r'$\kappa_\mathrm{L} \; \mathrm{(W/m/K)}$', fontsize=16)
##plt.plot(k_mags, kappa_list)
##plt.savefig('twistKappa_D1e-9.pdf', dpi=400, bbox_inches = 'tight')
##plt.show(block=False)
##    
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1}}$')
#plt.ylabel(r'$TBC$')
#plt.plot(k_mags, tbc_list)
#plt.savefig('tiltBoundary_tbc.pdf', dpi=400, bbox_inches = 'tight')
#plt.show()


'''
kappa versus T curve
'''
kappaT = []
kappaT_strain = []
kapitzaT = []
rk_T = []
temps = np.linspace(100, 300, 3)
i = 0
for T in temps:
    kappaT.append(AS.kL_T(Gamma_rot, k_max, dk, vg_k, omega_k, T, 50))
    rk_T.append(1/AS.tbc_T(k_max, dk, vg_k, omega_k, T, n_1D, Gamma_rot, 50))
#    kappaT_strain.append(AS.kL_T(Gamma, k_max, dk, vg_k, omega_k, T, 50))
#    kapitzaT.append(rk_T[i]*kappaT[i])
    i = i+1
plt.figure()
plt.plot(temps, kappaT)



'''
Save the results from the twist boundary calculation
'''
outfile = str(twist) + '_twist_vsT.npz'

np.savez(outfile, np.array(temps), np.array(kappaT), np.array(rk_T))

#%%
'''
Plots of thermal boundary conductance versus T
'''
#Thermal Bopundary Resistance figure
plt.figure()
plt.plot(temps, rk_T)
plt.xlabel('T (K)')
plt.ylabel(r'$R_K$ $(m^2K/W)$')
plt.savefig('tiltBoundary_Rk_T.pdf', bbox_inches = 'tight')

