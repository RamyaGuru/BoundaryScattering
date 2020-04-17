#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:09:04 2020

@author: ramyagurunathan

Helper functions
"""

from ThermalModel import Elasticity, CrystStruct, ThermalTransport
from math import pi as pi

'''
Constants
'''
hbar = 6.6260696 * 10**-34 / (2 * pi)
kB = 1.38064852e-23

def input_dict_from_MP(mpid):
    elast = Elasticity(mpid)
    struct = CrystStruct(mpid)
    thermal = ThermalTransport(struct, elast)
    input_dict = {'atmV' : struct.vol / struct.natoms,
                   'N' : struct.natoms,
                   'avg_vs' : thermal.therm_dict['average v_s'],
                   'gruneisen' : thermal.therm_dict['gruneisen'],
                   'nu' : elast.poisson,
                   }
    return input_dict
