from ArrayScattering import ArrayScattering as AS
import ScatteringPlots as SPlt
import ThermalTransport as TT
import math
import numpy as np
import helper
import json



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
#tilt = AS.ArrayScattering(**input_dict, geom = 'tilt', theta = 10, ax = 2, d_GS = 350E-9)

#Instantiate as object of ArrayScattering
tilt = AS(**input_dict, geom = 'tilt', theta = 10, ax = 2, d_GS = 350E-9)

# Scattering matrix elements
'''
q_vector[0] and q_vector[1] 
'''
def V1_twiddle_sq_Delta(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               ((tilt.b * (1 - 2 * tilt.nu)) / (1 - tilt.nu)) * (q_vector[1]\
               / (q_vector[0]**2 + q_vector[1]**2)))**2
# missing a negative sign from the factor of "i"?
    

def V1_twiddle_sq_S(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               (tilt.b / (1 - tilt.nu)) * ((q_vector[0] * q_vector[1]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2


def V1_twiddle_sq_R(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(helper.hbar * tilt.omega_kmag(k) * tilt.gruneisen * \
               tilt.b * ((2 * q_vector[0]) / (q_vector[0]**2 + q_vector[1]**2)))**2
#for the twist boundary case, the gruneisen parameter is unsed in the rotation term? Why?

def Gamma_GBS(k_vector, kprime_vectors):
   return tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta) \
          + tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S)\
          + tilt.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_R)

#Move to thermalTransport?
def Gamma(k_vector):
    return Gamma_GBS(k_vector, tilt.kprimes_y(k_vector)) * 1E-9 #what's this function for?

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

#Can this go in the ThermalTransport script?

def calculate_spectral_phonon_props(prop_list = ['tau', 'transmissivity', 'TBC', 'kappa'], n_angle = 100, n_k = 100, T = 300):
    '''
    Calculate frequency-dependent properties
    prop_list : spectral properties which should be calculated
    '''
    spectral = {'vg' : [], 'omega' : []}
    function = {'tau' : TT.tau_spectral, 'transmissivity' : TT.transmissivity_spectral,
                'TBC' : TT.tbc_spectral, 'kappa' : TT.kL_spectral}
    if any(prop_list) not in ['tau', 'transmissivity', 'TBC', 'kappa']:
        ValueError('Property not in allowed values list')
    if 'tau' in prop_list:
        spectral['tau'] = []
    if 'transmissivity' in prop_list:
        spectral['transmissivity'] = []
    if 'TBC' in prop_list:
        spectral['TBC'] = []
    if 'kappa' in prop_list:
        spectral['kappa'] = []
    dk = tilt.k_max / n_k
    k_mags = np.arange(dk, tilt.k_max, dk)
    params = {'Gamma' : Gamma, 'gb' : tilt, 'k' : dk, 'n_angle' : n_angle, 'T' : T}
    for k in k_mags: 
        spectral['vg'].append(tilt.vg_kmag(k))
        spectral['omega'].append(tilt.omega_kmag(k))
        params['k'] = k
        for prop in prop_list:
            spectral[prop].append(function[prop](**params))
    return spectral

if __name__ == "__main__":
    Gamma_list = calculate_Gammas(200)
    SPlt.diffraction_plot(tilt, Gamma_list[0], Gamma_list[1])
    SPlt.convergence_tau_plot(tilt, Gamma, 100, T = 300)
#    spectral = calculate_spectral_phonon_props(prop_list = ['tau', 'transmissivity', 'TBC', 'kappa'],\
#                                         n_angle = 100, n_k = 100, T = 300)
    with open('spectral.json') as json_file:
        spectral = json.load(json_file)
    SPlt.spectral_plots(tilt, spectral)



#%%Temperature plots

#%%
#rk_T = []
#kappaT = []
#kapitzaT = []
#temps = np.linspace(100, 500, 10)
#i=0
#for T in temps:
#    rk_T.append(1/AS.tbc_T(k_max, dk, vg_k, omega_k, T, n_1D, Gamma, 50))
#    kappaT.append(1/AS.kL_T(Gamma, k_max, dk, vg_k, omega_k, T, 50))
#    kapitzaT.append(rk_T[i]*kappaT[i])
#    i = i+1
##Thermal Bopundary Resistance figure
#plt.figure()
#plt.plot(temps, rk_T)
#plt.xlabel('T (K)')
#plt.ylabel(r'$R_K$ $(m^2K/W)$')
#plt.savefig('tiltBoundary_Rk_T.pdf', bbox_inches = 'tight')
#
##Thermal Conductivity figure
#plt.figure()
#plt.plot(temps, kappaT)
#plt.xlabel('T (K)')
#plt.ylabel(r'$\kappa_\mathrm{L} \; \mathrm{(W/m/K)}$')
#plt.savefig('tiltBoundary_kappa_T.pdf', bbox_inches = 'tight')
##Kapitza length figure
##plt.figure()
##plt.plot(temps, kapitzaT)
##plt.xlabel('T (K)')
##plt.ylabel(r'$L_K$ $m$')
##plt.savefig('tiltBoundary_kapitza_T.pdf', bbox_inches = 'tight')
###Log plots together with Qing Hao's Data
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