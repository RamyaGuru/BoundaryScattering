#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:40:45 2019

@author: ramyagurunathan

Acoustic Mismatch Model with the Elastic Anisotropy Term
"""

from pymatgen.ext.matproj import MPRester
import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/christoffel')
from christoffel import Christoffel
import numpy as np

mpr = MPRester('pfSJBa1OwitR5uNL')


'''
Structure: Crystal Structure
'''

class CrystStruct:
    def __init__(self, mp_id):
        self.struct = mpr.get_structure_by_material_id(mp_id)
        self.vol = self.struct.volume
        self.natoms = self.struct.composition.num_atoms
        self.nsites = self.struct.num_sites
        self.elem = self.struct.species
        self.speciesMass = self.struct.composition.weight
        self.atmMass = self.speciesMass/self.natoms

    def massDensity(self):
        return 1.6605E3 * self.nsites * self.speciesMass /\
    (self.natoms * self.vol)


'''
Material Property: Elasticity
'''

class Elasticity:
    def __init__(self, mp_id):
        self.xtal = CrystStruct(mp_id)
        self.elast = mpr.get_data(mp_id)
        self.bulkMod = self.elast[0]['elasticity']['K_VRH'] #Bulk modulus
        self.shearMod = self.elast[0]['elasticity']['G_VRH'] #shear modulus
        self.aniso = self.elast[0]['elasticity']['elastic_anisotropy'] #universal anisotropy value        
        self.stiffness = np.array(self.elast[0]['elasticity']['elastic_tensor']) #stiffness tensor
        self.density = self.xtal.massDensity()
        self.vp = {}
        self.vg = {}
        
    def get_group_phase_velocities(self, direction):
        '''
        Given the Materials Project elasticity and a direction, calculate the phase velocity
        '''
        chObj = Christoffel(self.stiffness, self.density)
        chObj.set_direction_cartesian(direction)
        dirkey = ''.join([n for n in direction])
        self.vp[dirkey] = chObj.get_phase_velocity()
        self.vg[dirkey] = chObj.get_group_velocity()
    
    def get_delta_q_tilt(theta):
        '''
        Convert the misorientation angle to a change in the direction
        '''
     
    def get_delta_q_twist():
        '''
        Convert misorientation angle to a change in direction
        '''
      
                
mp_id = 'mp-149'
Si_el = Elasticity(mp_id)

Si_el.get_group_phase_velocities([1,0,0])

