'''
General Array Scattering Methods

These are for an interface perpendicular to the x-axis
'''

import math
import numpy as np
import matplotlib.pyplot as plt
from helper import * 

'''
Constants
'''
hbar = 6.6260696 * 10**-34 / (2 * math.pi)
kB = 1.38064852e-23



class ArrayScattering:
    '''   
    Inputs:
       theta: GB angle in degrees
       vs: avg. sound velocity
       V: avg. atomic volume
       N: number of atoms per unit cell
       d_GS: avg. grain size
       ax: axis of dislocation spacing, x=0, y=1, z=2, dictionary for 'twist' or heterointerface'
       single value for 'tilt', none for twin
       geom = ['twist', 'tilt', 'twin']
       Initialize an ArrayScattering object with the crystal and microstructure inputs
    '''
    def __init__(self, avg_vs, atmV: list, N, d_GS, nu, theta = None, ax = None, geom = 'tilt', amm = None, gruneisen = 1):
       if geom not in ['twist', 'tilt', 'heterointerface', 'twin']:
           raise ValueError('GB geometry value not valid')
       self.theta = theta
       self.vs = avg_vs
       self.V = atmV[0]
       self.N = N
       self.d_GS = d_GS
       self.n_1D = 3/d_GS
       self.gruneisen = gruneisen
       self.b = (self.V * self.N) ** (1 / 3) 
       self.k_max = (6 * math.pi**2 / (self.V * self.N))**(1 / 3)
       self.nu = nu
       self.omegaD = self.vs * self.k_max
       self.geom = geom
       if amm:
           self.amm = amm
       if geom == 'twist' or geom == 'tilt':
           self.theta = theta
           self.D = self.b / ( 2 * math.tan(theta * (math.pi/180) / 2))
           self.ax = ax 
       elif geom == 'twin':
           self.theta = theta
       elif geom == 'heterointerface':
           self.V2 = atmV[1]
           self.D = (self.V**(1/3) / (self.V2**(1/3) - self.V**(1/3))) * self.b
           self.ax = ax
                  
    
    '''
    Dispersion Properties:
        Need to make these more general
    '''
    def k_mag(k_vector):
        '''
        Input: k_vector expressed as array or list of length 3
        Output: scalar magnitude value
        Description:
        Outputs magnitude of the k_vector
        '''
        return (k_vector[0]**2 + k_vector[1]**2 + k_vector[2]**2)**(1 / 2)
    
    def gamma1(grun):
        '''
        Input: k_vector
        Output: gruneisen parameter as a scalar value
        '''
        return grun
    
    def omega_kmag(self, kmag):
        '''
        Input: k_vector
        Output: freqeuncy (scalar)
        '''
        return self.vs * kmag
    
    def vg_kmag(self, kmag):
        '''
        Input: k_vector
        Output: group velocity (scalar)
        '''  
        return self.vs
    
    '''
    Read Shockley Grain Boundary Energy
    '''
#    def gb_energy(self):
        
    
    '''
    Scattering Functions
    ''' 
    def qm(self, m):
        '''
        Input: index m and dislocation spacing D
        Output: scattering vector q_m
        Description:
        Scattering vector of index m
        '''
        return 2 * math.pi * m / self.D
    
    
    def kxprime_msigma(self, kx, kD, m, sign):
        '''
        Input: 
            kx: initial x-component of k-vector
            kD: component of k-vector along axis of dislocation spacing
            m: index of scattering vector
            sign: one of the two different final vectors (+/-)
            D: dislocation spacing
        Output: kx' for a specific scattering vector q_m
        kxprime derived by asserting the conservation laws:
        k'=k, kz'=kz, and ky'=ky-qm.
        '''
        return (sign * (kx**2 + (2 * kD * self.qm(m)) - self.qm(m)**2)**(1 / 2))
    
    
    def qx_msigma(self, kx, ky, m, sign):
        '''
        Input: 
            kx, ky: x and y component of the k-vector
            m: index of the scattering vector
            sign: pick out +/- final vector
            D: dislocation spacing
        Output: x-component of the scattering vector: qx = kxprime - kx
        '''
        return self.kxprime_msigma(kx, ky, m, sign, self.D) - kx
    
    
    def kprimes_y(self, k_vector):
        '''
        Input:
            k_vector
            D: dislocation spacing
        Output: list of k'_y values
        Calculates all possible kprime_vectors for a given k_vector.
        This applies the conservation laws, k'=k, kz'=kz, and ky'=ky-qm.
        '''
        kx = k_vector[0]
        ky = k_vector[1]
        kz = k_vector[2]
    
        # largest poisive m value set by requiring that
        # (kx**2 + (2 * ky * qm(m, D)) - qm(m, D)**2) >= 0.
        m_maxPos = math.modf((self.D * (ky + (kx**2 + ky**2)**(1 / 2))) / (2 * math.pi))[1]
        # largest negative m value set by the same condition
        m_maxNeg = math.modf((self.D * (ky - (kx**2 + ky**2)**(1 / 2))) / (2 * math.pi))[1]
    
        # make list of possible m values which are non-zero... so this is ignoring m=0 term?
        m_values = []
        for m in range(int(m_maxNeg), 0):
            m_values.append(m)
        for m in range(1, int(m_maxPos) + 1):
            m_values.append(m)
    
        # make list of possible k' vectors
        kprime_list = []
        # Add back-scattering, sign=-1 and m=0. Omit forward scattering
        # sign = 1 and m=0 or (k_vector=kprime_vector) as this causes a
        # 0/0 indeterminate in the V1_twiddle_sqs.
        kprime_list.append([-kx, ky, kz])
        for m in m_values:
            for sign in [-1, 1]:
                kprime_list.append([self.kxprime_msigma(kx, ky, m, sign), ky - self.qm(m), kz])
    
        return kprime_list
    
    def kprimes_z(self, k_vector):
        '''
        Input:
            k_vector
            D: dislocation spacing
        Output: list of k'_y values
        Calculates all possible kprime_vectors for a given k_vector.
        This applies the conservation laws, k'=k, kz'=kz, and ky'=ky-qm.
        '''
        kx = k_vector[0]
        ky = k_vector[1]
        kz = k_vector[2]
    
        # largest poisive m value set by requiring that
        # (kx**2 + (2 * ky * qm(m, D)) - qm(m, D)**2) >= 0.
        m_maxPos = math.modf((self.D * (kz + (kx**2 + kz**2)**(1 / 2))) / (2 * math.pi))[1]
        # largest negative m value set by the same condition
        m_maxNeg = math.modf((self.D * (kz - (kx**2 + kz**2)**(1 / 2))) / (2 * math.pi))[1]
    
        # make list of possible m values which are non-zero
        m_values = []
        for m in range(int(m_maxNeg), 0):
            m_values.append(m)
        for m in range(1, int(m_maxPos) + 1):
            m_values.append(m)
    
        # make list of possible k' vectors
        kprime_list = []
        # Add back-scattering, sign=-1 and m=0. Omit forward scattering
        # sign = 1 and m=0 or (k_vector=kprime_vector) as this causes a
        # 0/0 indeterminate in the V1_twiddle_sqs.
        kprime_list.append([-kx, ky, kz])
        for m in m_values:
            for sign in [-1, 1]:
                kprime_list.append([self.kxprime_msigma(kx, kz, m, sign), ky, kz-self.qm(m)])
    
        return kprime_list
    
    def kprimes_plot(self, k_vector):
        '''
        Input:
        Output: 
        Plots the scattering diagram. (WORKS BUT NEEDS SOME LOVE TO MAKE IT PRETTY)
        '''
        k = ArrayScattering.k_mag(k_vector)
    
        dtheta = math.pi / 100.
        kx_cirlist = []
        ky_cirlist = []
        for theta in np.arange(0, 2 * math.pi + dtheta, dtheta):
            kx_cirlist.append(k * math.cos(theta))
            ky_cirlist.append(k * math.sin(theta))
        plt.figure()
        plt.plot(kx_cirlist, ky_cirlist)
        plt.plot(k_vector[0], k_vector[1], 'ko', markersize=8)
        plt.plot(np.asarray(self.kprimes_y(k_vector))[:, 0], \
                 np.asarray(self.kprimes_y(k_vector))[:, 1], 'ro', markersize=4)
        plt.axes().set_aspect(aspect=1)
        plt.show(block=False)
    
    
    def GammaArray(self, k_vector, kprime_vectors, V1_twiddle_sq, ax):
        '''
        Performs sum over all possible k' states.
        Requires the magnitude squared of the Fourier transform
        of the scattering potential, V1Twidle2(k_vector, kprime_vector).
        ax: axis of the dislocation spacing
        
        Double-checked: 1/7/2020
        Rather than just Gamma, this is basically relaxation time?
        '''
        k = ArrayScattering.k_mag(k_vector)
        kx = k_vector[0]
        kD = k_vector[ax]
        # sum over all possible kprime_vectors
        i=0
        running_sum = 0
        for kprime_vector in kprime_vectors:
            kxprime = kprime_vector[0]
            kDprime = kprime_vector[ax]
            qx = kxprime - kx
            qD = kDprime - kD
            running_sum = running_sum + \
            V1_twiddle_sq(self, k_vector, kprime_vector) * (-qx * kx - qD * kD) * abs(kxprime) ** -1  
            i+=1
        return (self.n_1D / (hbar ** 2 * self.D ** 2 * self.vg_kmag(k) * k)) * running_sum  
    
    def GammaArray_rot(self, k_vector, V_twiddle_R):
        k = ArrayScattering.k_mag(k_vector)
        gamma = (self.n_1D / (hbar ** 2 * self.vs)) * V_twiddle_R(self, k_vector)
        return gamma
    