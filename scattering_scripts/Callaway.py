import math
import numpy as np

kB = 1.3806485E-23
hbar = 6.6260696 * 10**-34 / (2 * math.pi)


def cv_ph(omega, T):
    x = (hbar * omega) / (kB * T)
    return (kB * x**2 * np.exp(x)) / (np.exp(x) - 1)**2


def kappa_L_discrete(k_list, omega_list, vg_list, tau_list, T):
    running_integrand = 0
    dk = k_list[1] - k_list[0]
    i = 0
    for k in k_list:
        omega = omega_list[i]
        vg = vg_list[i]
        tau = tau_list[i]
        running_integrand = running_integrand + cv_ph(omega, T) * vg**2 * tau * k**2
    return (1 / (2 * math.pi**2)) * running_integrand * dk
