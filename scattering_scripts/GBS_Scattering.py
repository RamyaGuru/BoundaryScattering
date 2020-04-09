import ArrayScattering as AS
import Callaway as Cal
import math
import numpy as np
import matplotlib.pyplot as plt


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
    k = AS.k_mag(k_vector)
    return vs * k
def vg_k(k_vector):
    return vs


# Microstructural parameters
b = (V * N) ** (1 / 3)  # Burger's vector [m]
d_GS = 350E-9           # Average grain size [m]
n_1D = 3 / d_GS         # Number density of GBs [m^-1]
D = 1E-9                # Linear defects spacing [m]


# Scattering matrix elements
'''
q_vector[0] and q_vector[1] 
'''
def V1_twiddle_sq_Delta(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(AS.hbar * omega_k(k_vector) * gamma(k_vector) * \
               ((b * (1 - 2 * nu)) / (1 - nu)) * (q_vector[1]\
               / (q_vector[0]**2 + q_vector[1]**2)))**2
# missing a negative sign from the factor of "i"?
    

def V1_twiddle_sq_S(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(AS.hbar * omega_k(k_vector) * gamma(k_vector) * \
               (b / (1 - nu)) * ((q_vector[0] * q_vector[1]**2)\
               / (q_vector[0]**2 + q_vector[1]**2)**2)) ** 2


def V1_twiddle_sq_R(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(AS.hbar * omega_k(k_vector) * gamma(k_vector) * \
               b * ((2 * q_vector[0]) / (q_vector[0]**2 + q_vector[1]**2)))**2
#for the twist boundary case, the gruneisen parameter is unsed in the rotation term? Why?

def Gamma_GBS(k_vector, kprime_vectors, vg, n_1D, D):
   return AS.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_Delta, vg, n_1D, D, 1) \
          + AS.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_S, vg, n_1D, D, 1)\
          + AS.GammaArray(k_vector, kprime_vectors, V1_twiddle_sq_R, vg, n_1D, D, 1)


def Gamma(k_vector, vg):
    return Gamma_GBS(k_vector, AS.kprimes_y(k_vector, D), vg, n_1D, D) * 1E-9 #what's this function for?

k_vector = [0.2 * k_max, 0, 0]



# Plot Gamma_GBS(k_vector) for normal incidence
n_k = 20
dk = k_max / n_k
k_mags = np.arange(dk, k_max, dk)
k_norm = k_mags / k_max
k_vectors = []
for k in k_mags:
    k_vectors.append([k, 0, 0])

Gamma_GBS_list = []
for k_vector in k_vectors:
    Gamma_GBS_list.append(Gamma(k_vector, vg_k(k_vector)))

plt.figure()
plt.xlim((0, 1))
plt.ylim((0, 40))
plt.xlabel(r'$k/k_{\mathrm{max}}$', fontsize=16)
plt.ylabel(r'$\Gamma \; \mathrm{(ns^{-1})}$', fontsize=16)
plt.plot(k_norm, Gamma_GBS_list)
plt.savefig('tiltdiff_D1e-9_2.pdf', dpi=400, bbox_inches = 'tight')
plt.show(block=False)

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
omega_list = []
vg_list = []
tau_list = []
trans_list = []
tbc_list = []
kappa_list = []

#%% Spectral plots
#NOTE: already assume that all the phonons are incident to the interface
T = 300
for k in k_mags:
    omega_list.append(omega_k([k,0,0])) # omega and vg are supposed to be only a function of k, not k_vector. This is tacky and needs to be fixed!
    vg_list.append(vg_k([k,0,0]))
    tau_list.append(AS.tau_spectral(Gamma, k, vg_k, 50))
    trans_list.append(AS.transmissivity(k, vg_k, n_1D, Gamma, 50))
    tbc_list.append(AS.tbc_spectral(k, vg_k, omega_k, T, Gamma, n_1D, 50))
    kappa_list.append(AS.kL_spectral(Gamma, k, vg_k, omega_k, T, 50))
#
#
plt.figure()
plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)
plt.plot(k_mags, tau_list)
plt.savefig('tiltBoundary_D1e-9_2.pdf', dpi=400, bbox_inches = 'tight')
plt.show(block=False)
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