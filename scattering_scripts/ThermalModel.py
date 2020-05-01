# -*- coding: utf-8 -*-
"""
Spyder Editor

Script to work with Materials Project REST API
Try out a Boltzmann Transport Equation calc
"""

from pymatgen.ext.matproj import MPRester

mpr = MPRester('pfSJBa1OwitR5uNL')

#bands = mpr.get_bandstructure_by_material_id("mp-19717", line_mode = False)


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
        self.elast = mpr.get_data(mp_id) 
        self.bulkMod = self.elast[0]['elasticity']['K_VRH'] #Bulk modulus
        self.shearMod = self.elast[0]['elasticity']['G_VRH'] #shear modulus
        self.poisson = self.elast[0]['elasticity']['homogeneous_poisson'] # poisson ratio
        self.elastic_tensor = self.elast[0]['elasticity']['elastic_tensor'] # elastic tensor
        if self.elast[0]['elasticity']['warnings'] is not None:
            print(self.elast[0]['elasticity']['warnings'])
            
'''
Transport Property: Thermal Conductivity Calculation
'''

class ThermalTransport:
    #Longitudinal speed of sound calculation
    def __init__(self, structure, elasticity):
        self.matstruct = structure
        self.matprops = elasticity
        self.temperature = 300
        self.therm_coeffs = [6e-8, 0.012]
        self.therm_dict = {
            'bulk modulus' : self.matprops.bulkMod,
            'shear modulus' : self.matprops.shearMod,
            'average v_s' : self.avg_vs(),
            'gruneisen' : self.grunProxy2(),
            'kappa_lat' : self.kappa_lat(300, 12e-9, 0.12),
        }
    
    #Longitudinal speed of sound
    def trans_v(self):
        density = self.matstruct.massDensity()
        return (1E9 * self.matprops.shearMod / density)**0.5
    #Transverse speed of sound
    def long_v(self):
        density = self.matstruct.massDensity()
        return (1E9 * (self.matprops.bulkMod + 4./3. * self.matprops.shearMod) /\
                density)**0.5
    #Speed of sound calculation
    def avg_vs(self):
        vl = self.long_v()
        vt = self.trans_v()
        return 3**(1/3)*(1/vl**3 + 2/vt**3)**(-1/3)
    

    def grunProxy(self):
        vl = self.long_v()
        vt = self.trans_v()
        return vt / vl

    #Gruneisen approximation: full expression from M. Agne
    def grunProxy2(self):
        vl = self.long_v()
        vt = self.trans_v()
        gamma = (9/2)*((vl/vt)**2 - (4/3))/((vl/vt)**2 + 2)
        return gamma
    #Acoustic thermal conductivity contribution
    def kappa_a(self, T):
        vs= self.avg_vs()
        grun = self.grunProxy2()
        s = self.matstruct
        return float(s.atmMass) * vs**3 / (T * s.vol**(2/3) * grun**2) *\
    (1 / s.natoms**(1/3))
    
    #Optical thermal conductivity contribution
    def kappa_o(self):
        vs = self.avg_vs()
        return (vs/ self.matstruct.vol**(2/3)) *\
    (1 - self.matstruct.natoms**(-2/3))

    #Need to perform regression to get the rest of the properties
    '''
    kappa_tot method that combines kappa_a and kappa_o with the fitting
    parameters
    '''
    def kappa_lat(self, T, A, B):
        return A*self.kappa_a(T) + B*self.kappa_o()
'''
Transport Transport Property: Point Defect Scattering Effects on Thermal Conductivity
'''    

     
