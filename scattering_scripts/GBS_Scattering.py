from ArrayScattering import ArrayScattering as AS
import ScatteringPlots as splt
import Callaway as Cal
import math
import numpy as np
import matplotlib.pyplot as plt
import helper



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
    return [list(a) for a in zip(Gamma_GBS_list, k_norm)]

Gamma_GBS_list = calculate_Gammas(20)


# Plot Gamma_GBS(k_vector) for normal incidence
#n_k = 20
#dk = k_max / n_k
#k_mags = np.arange(dk, k_max, dk)
#k_norm = k_mags / k_max
#k_vectors = []
#for k in k_mags:
#    k_vectors.append([k, 0, 0])
#
#Gamma_GBS_list = []
#for k_vector in k_vectors:
#    Gamma_GBS_list.append(Gamma(k_vector, vg_k(k_vector)))

#plt.figure()
#plt.xlim((0, 1))
#plt.ylim((0, 40))
#plt.xlabel(r'$k/k_{\mathrm{max}}$', fontsize=16)
#plt.ylabel(r'$\Gamma \; \mathrm{(ns^{-1})}$', fontsize=16)
#plt.plot(k_norm, Gamma_GBS_list)
#plt.savefig('tiltdiff_D1e-9_2.pdf', dpi=400, bbox_inches = 'tight')
#plt.show(block=False)

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



# Calculation of spectral tau and kappa
#omega_list = []
#vg_list = []
#tau_list = []
#trans_list = []
#tbc_list = []
#kappa_list = []
#
##%% Spectral plots
##NOTE: already assume that all the phonons are incident to the interface
#T = 300
#for k in k_mags:
#    omega_list.append(omega_k([k,0,0])) # omega and vg are supposed to be only a function of k, not k_vector. This is tacky and needs to be fixed!
#    vg_list.append(vg_k([k,0,0]))
#    tau_list.append(AS.tau_spectral(Gamma, k, vg_k, 50))
#    trans_list.append(AS.transmissivity(k, vg_k, n_1D, Gamma, 50))
#    tbc_list.append(AS.tbc_spectral(k, vg_k, omega_k, T, Gamma, n_1D, 50))
#    kappa_list.append(AS.kL_spectral(Gamma, k, vg_k, omega_k, T, 50))
##
##
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
#plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)
#plt.plot(k_mags, tau_list)
#plt.savefig('tiltBoundary_D1e-9_2.pdf', dpi=400, bbox_inches = 'tight')
#plt.show(block=False)
#
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1}}$')
#plt.ylabel(r'$t$')
#plt.plot(k_mags, trans_list)
#plt.savefig('tiltBoundary_trans.pdf', dpi=400, bbox_inches = 'tight')
#plt.show()
#
#
#plt.figure()
#plt.xlabel(r'$k \; \mathrm{(m^{-1}}$')
#plt.ylabel(r'$TBC$')
#plt.plot(k_mags, tbc_list)
#plt.savefig('tiltBoundary_tbc.pdf', dpi=400, bbox_inches = 'tight')
#plt.show()

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