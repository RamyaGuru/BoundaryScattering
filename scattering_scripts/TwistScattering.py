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
    
    need to figure out ax for twist boundary script

These are (001) type twist grain boundaries
"""

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import ScatteringPlots as SPlt
from AngularVs import ElasticProps, v_long, v_shear, v_sound, v_xy
import math
import numpy as np
from math import asin, acos, pi, sin, cos, tan, atan
import helper
from AMMTransport import AMMTransport

np.seterr(divide='raise', invalid="raise")


cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330


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


    

'''
Input dict
'''
input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

'''
Initialize input dictionary with Materials Project
'''
#input_dict = helper.input_dict_from_MP('mp-149')

twist = AS(**input_dict, geom = 'twist', theta = 5, ax = {'n' : 1, 'm' : 2}, d_GS = 350E-9)
amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])
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

'''
Alternate rotation potentials (previous tries)
'''    

def V_twiddle_sq_R(k_vector):
#    k = AS.k_mag(k_vector)
#    q_vector = np.asarray(kprime_vector)- np.asarray(k_vector)
#    #calculate the angle of incidence, what about misorientation angle?
#    inc = acos(k_vector[0]/(k))
    vs_new = v_sound(twist.theta * (math.pi / 180))
    # k will cancel out
    return (helper.hbar*(abs(vs_new - twist.vs)/(2*math.pi)))**2  #multiply by the ratio of k to q?

#def V_R_ch(k_vector, kprime_vector):
#    k = AS.k_mag(k_vector)
#    knorm = k_vector / k
#    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
#    v1, v2 = amm.vs_rot(knorm, helper.rot_tensor_x, twist.theta)
#    return ((helper.hbar / (2 * pi)) * abs(1000 * (v2 - v1)) * (k / q_vector[0]))**2  
    

#Shouldn't the rotation be about x?
def V_R_ch_snell(k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
#    theta1 = acos(k_vector[0]/kmag)
    v1, v2, theta2 = amm.vs_rot_Snell(knorm, helper.rot_tensor_x, twist.theta)
#    print([v1, v2])
#    if math.isnan(theta2):
#        qx = -2 * k_vector[0]
#        q_vector = [qx, 0, 0]
#    else:
#        qx = kmag * cos(theta2) - k_vector[0]
#        q_vector = [qx, (k_vector[1] / sin(theta1)) * sin(theta2) - k_vector[1], (k_vector[2] / sin(theta1)) * sin(theta2) - k_vector[2]]
#    print(-1 * np.dot(q_vector, k_vector))
#    print((v2 - v1) / v1)
    return abs(helper.hbar * (abs(v2 - v1) / v1) * twist.vs * (kmag / (2 * k_vector[0])))**2 * (2 * k_vector[0]**2 / kmag**2)
#    return abs(helper.hbar * (abs(v2 - v1) / v1) * twist.vs)**2 * (2 * k_vector[0]**2 / kmag**2)
#    return abs(helper.hbar * (abs(v2 - v1) / v1) * (twist.omega_kmag(kmag) / (2 * k_vector[0])))**2 * (2 * k_vector[0]**2 / kmag**2) #replaced qx by 1/ Burger's vector (width of delta function scattering potential)


'''
Rotation Potential
'''

def V_tilde_sq_R(k_vector):
    kmag = helper.k_mag(k_vector)
    knorm = k_vector / kmag
    v1, v2, theta2 = amm.vs_rot_Snell(knorm, helper.rot_tensor_x, twist.theta)
    return abs(helper.hbar * (abs(v2 - v1) / v1) * twist.vs * (kmag / (2 * k_vector[0])))**2 * (2 * k_vector[0]**2 / kmag**2)

'''
Strain components of the "n" array: Dislocaiton line in z, spacing in y
'''

def Vn_twiddle_sq_E13(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * twist.omega_kmag(k)* twist.gruneisen *(twist.b*q_vector[1])/(2* (q_vector[0]**2 + q_vector[1]**2)))**2
    
def Vn_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * twist.omega_kmag(k) * twist.gruneisen * (twist.b * q_vector[0]/(2*(q_vector[0]**2 + q_vector[1]**2))))**2

def V_twiddle_sq_n(k_vector, kprime_vector):
    return Vn_twiddle_sq_E13(k_vector, kprime_vector) + Vn_twiddle_sq_E23(k_vector, kprime_vector)
    

'''
Strain components of the "m" array: Dislocation line in y, spacing in z
'''
def Vm_twiddle_sq_E12(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * twist.omega_kmag(k) * twist.gruneisen * (twist.b*q_vector[2])/(2* (q_vector[0]**2 + q_vector[2]**2)))**2
    
def Vm_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * twist.omega_kmag(k) * twist.gruneisen *(twist.b*q_vector[0]/(2*(q_vector[0]**2 + q_vector[2]**2))))**2

def V_twiddle_sq_m(k_vector, kprime_vector):
    return Vm_twiddle_sq_E12(k_vector, kprime_vector) + Vm_twiddle_sq_E23(k_vector, kprime_vector)

'''
STGB scattering rates
'''
#Still need to include the rotation term
def Gamma_GBS(k_vector, kprime_yvectors, kprime_zvectors):
   return twist.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, twist.ax['n']) \
          + twist.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, twist.ax['m']) 


def Gamma_GBS_rot(k_vector, kprime_yvectors, kprime_zvectors):
   tot = [twist.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, twist.ax['n']) \
          ,twist.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, twist.ax['m']) \
          ,twist.GammaArray_rot(k_vector, V_tilde_sq_R)]
   return sum(tot)

def Gamma_GBS_rot_only(k_vector, kprime_yvectors, kprime_zvectors):
    return twist.GammaArray_rot(k_vector, V_tilde_sq_R) #twist.GammaArray_rot(k_vector, V_R_ch_snell)
          
def Gamma(k_vector):
    return Gamma_GBS(k_vector, twist.kprimes_y(k_vector), twist.kprimes_z(k_vector)) * 1E-9 

def Gamma_rot(k_vector):
    return Gamma_GBS_rot(k_vector, twist.kprimes_y(k_vector), twist.kprimes_z(k_vector)) * 1E-9

def Gamma_rot_only(k_vector):
    return Gamma_GBS_rot_only(k_vector, twist.kprimes_y(k_vector), twist.kprimes_z(k_vector)) * 1E-9

def calculate_Gammas(n_k):
    dk = twist.k_max / n_k
    k_mags = np.arange(dk, twist.k_max, dk)
    k_norm = k_mags / twist.k_max
    k_vectors = []
    for k in k_mags:
        k_vectors.append([k, 0, 0]) #all same direction...
    Gamma_GBS_list = []
    for k_vector in k_vectors:
        Gamma_GBS_list.append(Gamma_rot(k_vector))
    return [k_norm, Gamma_GBS_list]    


if __name__ == "__main__":
#    Gamma_list = calculate_Gammas(200)
#    SPlt.diffraction_plot(twist, Gamma_list[0], Gamma_list[1])
    SPlt.convergence_tau_plot(twist, Gamma_rot_only, 600, T = 300, save = True)
#    spectral = TT.calculate_spectral_props(twist, Gamma_rot_only, prop_list = ['tau'],\
#                                        n_angle = 200, n_k = 10, T = 300) #n_angle = 200, n_k = 100
#    SPlt.spectral_plots(twist, spectral, prop_list = ['tau'], save = True)
#    temp_dependence = TT.calculate_temperature_dependence(twist, Gamma, temp_list = [100, 800])



'''
kappa versus T curve
'''
#kappaT = []
#kappaT_strain = []
#kapitzaT = []
#rk_T = []
#temps = np.linspace(100, 300, 3)
#i = 0
#for T in temps:
#    kappaT.append(AS.kL_T(Gamma_rot, k_max, dk, vg_k, omega_k, T, 50))
#    rk_T.append(1/AS.tbc_T(k_max, dk, vg_k, omega_k, T, n_1D, Gamma_rot, 50))
##    kappaT_strain.append(AS.kL_T(Gamma, k_max, dk, vg_k, omega_k, T, 50))
##    kapitzaT.append(rk_T[i]*kappaT[i])
#    i = i+1
#plt.figure()
#plt.plot(temps, kappaT)



'''
Save the results from the twist boundary calculation
'''
#outfile = str(twist) + '_twist_vsT.npz'
#
#np.savez(outfile, np.array(temps), np.array(kappaT), np.array(rk_T))
#
##%%
#'''
#Plots of thermal boundary conductance versus T
#'''
##Thermal Bopundary Resistance figure
#plt.figure()
#plt.plot(temps, rk_T)
#plt.xlabel('T (K)')
#plt.ylabel(r'$R_K$ $(m^2K/W)$')
#plt.savefig('tiltBoundary_Rk_T.pdf', bbox_inches = 'tight')

