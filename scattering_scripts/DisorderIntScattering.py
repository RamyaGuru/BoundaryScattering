#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 17:53:24 2019

@author: ramyagurunathan

Disorder Scattering Model from the Beechem 2007 paper 

Virtual Crystal Approximation Diffuse Mismatch Model
"""

import ArrayScattering as AS
import AlloyScattering as AY

hbar = 6.6260696 * 10**-34 / (2 * math.pi)
kB = 1.38064852e-23



'''
Define interface properties
'''
vs_1 = 6084 #Si speed of sound
vs_2 = 5400 #Ge speed of sound 
D = 1.2e-9 #Interface thickness in nanometers

'''
Define sample properties
'''


def vs_vca(c: int):
    '''
    Sound velocity for the virtual crystal at the interface
    
    c: the fractional composition of each component
    '''
    vs_vca = vs_1*c + vs_2*(1-c)
    return vs_vca

def alpha(vsA, vsB):
    '''
    Transmissivity from A -> B expression based on the DMM.
    
    Factors in the spectral density of states cancel, so that you just get a 
    ratio of inverse square speeds of sounds
    
    Do we need the geometric factor based on the angle of incidence?
    Should just give a factor of 1/2
    '''
    alpha = (1/2)*(vsB**(-2) / (vsA**(-2) + vsB**(-2)))
    return alpha

'''
Computation of delta: determines likelihood of multiple scattering events in the interface region
'''

#def vca_mfp():
#    '''
#    Mean free path in virtual crystal (obtain from the thermal conductivity value)
#    '''
    
    
def transmissivity(vsA, vsB, c):
    '''
    The alphabet soup mess means: 
        Thermal Boundary Conductance_ Virtual Crystal Approximation_ Diffuse Mismatch Model
    '''
    vsInt = vs_vca(c)
    alpha1 = alpha(vsA, vsInt)
    alpha2 = alpha(vs_vca, vsB)
    alpha = alpha1 + alpha2
    return alpha

    

    
