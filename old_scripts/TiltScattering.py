from ArrayScattering import ArrayScattering as AS
from AMMTransport import AMMTransport
import ScatteringPlots as SPlt
import ThermalTransport as TT
import math
import numpy as np
import helper
import json

'''
Tilt Scattering Riley Scripts
'''


'''
Input dict
'''
input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330

geom = 'tilt'

'''
Initialize input dictionary with Materials Project
'''
#input_dict = helper.input_dict_from_MP('mp-149')

#Instantiate as object of ArrayScattering
tilt = AS(**input_dict, geom = 'tilt', theta = 5, ax = 2, d_GS = 350E-9)


'''
Initialize function
'''
def initialize(input_dict, cmat, density, theta, geom, ax = 1, d_GS = 350e-9):
    amm = AMMTransport(cmat, density, input_dict['atmV'][0], input_dict['N'])
    tilt = AS(**input_dict, geom = geom, amm = amm, theta = theta, ax = ax, d_GS = d_GS)
    return tilt

# Scattering matrix elements
'''
q_vector[0] and q_vector[1] 
'''
def V1_twiddle_sq_Delta(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               ((tilt.b * (1 - 2 * tilt.nu)) / (1 - tilt.nu)) * (q_vector[1]\
               / (q_vector[0]**2 + q_vector[1]**2)))**2
# missing a negative sign from the factor of "i"?
    

def V1_twiddle_sq_S(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               (tilt.b / (1 - tilt.nu)) * ((q_vector[0] * q_vector[1]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2


def V1_twiddle_sq_R(tilt, k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               tilt.b * ((2 * q_vector[0]) / (q_vector[0]**2 + q_vector[1]**2)))**2
#for the twist boundary case, the gruneisen parameter is unsed in the rotation term? Why?

def Gamma_GBS(tilt, k_vector, kprime_vectors):
    tot = [tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, tilt.ax) \
          ,tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, tilt.ax)\
          ,tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_R, tilt.ax)]
    return sum(tot)

#Move to thermalTransport?
def Gamma(tilt, k_vector):
    return Gamma_GBS(tilt, k_vector, tilt.kprimes_y(k_vector)) * 1E-9 #what's this function for?

k_vector = [0.2 * tilt.k_max, 0, 0]


def calculate_Gammas(n_k):
    dk = tilt.k_max / n_k
    k_mags = np.arange(dk, tilt.k_max, dk)
    k_norm = k_mags / tilt.k_max
    k_vectors = []
    for k in k_mags:
        k_vectors.append([k, 0, 0])
    Gamma_GBS_list = []
    for k_vector in k_vectors:
        Gamma_GBS_list.append(Gamma(k_vector))
    return [k_norm, Gamma_GBS_list]    




if __name__ == "__main__":
#    Gamma_list = calculate_Gammas(200)
#    SPlt.diffraction_plot(tilt, Gamma_list[0], Gamma_list[1])
    theta = 5
    tilt = initialize(input_dict, cmat, density, theta, geom = 'tilt', ax = 1, d_GS = 350e-9)
    SPlt.convergence_tau_plot(tilt, Gamma, 110, T = 300)
#    spectral = TT.calculate_spectral_props(tilt, Gamma, prop_list = ['tau'],\
#                                         n_angle = 200, n_k = 100, T = 300)
#    with open('spectral.json') as json_file:
#        spectral = json.load(json_file)
#    SPlt.spectral_plots(tilt, spectral, prop_list = ['tau'], save = True)
#    temp_dependence = TT.calculate_temperature_dependence(tilt, Gamma, temp_list = [100, 800])


    
#%%
'''
Log plots together with Qing Hao's Data
'''
#
#qh_data = np.loadtxt('Si_tbr.csv', delimiter = ',')
#qh_data[:,1] = qh_data[:,1]*1e-9
##Data
#plt.figure()
#plt.scatter(qh_data[:,0], qh_data[:,1], label = 'QH data')
#plt.xscale('log')
#plt.yscale('log')
#
#plt.loglog(qh_data[:,0], qh_data[0,1]*(qh_data[:,0]/qh_data[0,0])**(-1), label = r'T$^{-1}$')
#plt.loglog(qh_data[:,0], qh_data[0,1]*(qh_data[:,0]/qh_data[0,0])**(-1.75), label = r'T$^{-1.75}$')
#plt.legend()
#plt.xlabel('T (K)')
#plt.ylabel(r'$R_K$ $(m^2K/W)$')
#plt.ylim([1e-9,1e-6])
#plt.savefig('qhEdepedence.pdf', bbox_inches = 'tight')
#
#
##Model
#plt.figure()
#plt.loglog(temps, rk_T, label = 'Tilt Boundary model')
#plt.xlabel('T (K)')
#plt.ylabel(r'$R_K$ $(m^2K/W)$')
#
##Fit
#plt.loglog(temps, rk_T[0]*(temps/temps[0])**(-1), label = r'T$^{-1}$')
#plt.loglog(temps, rk_T[0]*(temps/temps[0])**(-1.4), label = r'T$^{-1.4}$')
#plt.legend()
#plt.xlabel('T (K)')
#plt.ylabel(r'$R_K$ $(m^2K/W)$')
#plt.yscale()
#plt.ylim([1e-9,1e-6])
#plt.savefig('modelTdependence.pdf', bbox_inches = 'tight')