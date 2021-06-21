#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:07:19 2021

@author: ramyagurunathan

Study temperature dependence of the Grid model
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/BoundaryScattering/scattering_scripts')

from ArrayScattering import ArrayScattering as AS
import ThermalTransport as TT
import TwistScattering as TS
import time
import ScatteringPlots as SPlt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
import matplotlib.ticker as ticker
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec

cmat = np.array([[165.6, 63.9, 63.9, 0, 0, 0],[63.9, 165.6, 63.9, 0, 0, 0],[63.9, 63.9, 165.6, 0, 0 ,0],\
     [0, 0, 0, 79.5, 0, 0],[0, 0, 0, 0, 79.5, 0],[0 ,0, 0, 0, 0, 79.5]])

density = 2330 / 1.0018

#Updated input dictionary based on thermal expansion

input_dict = {'avg_vs': 6094,
             'atmV': [2E-29 * 1.0018],
             'N': 2,
             'bulkmod' : 97.83E9,
             'shearmod' : 60E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = TS.initialize(input_dict, cmat, density, geom = 'twist', theta = 12, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)


start_time = time.time()
tau_list1 = SPlt.convergence_tau_plot(twist, TS.Gamma_rot, 150, T = 300, save = True)
print("--- %s seconds ---" % (time.time() - start_time))


#Old input dictionary
density = 2330 / 1.0018

input_dict = {'avg_vs': 6084,
             'atmV': [2E-29],
             'N': 2,
             'bulkmod' : 97.83E9,
             'shearmod' : 60E9,
             'nu' : 0.29,
             'gruneisen' : 1,
        }

twist = TS.initialize(input_dict, cmat, density, geom = 'twist', theta = 12, ax = {'n': 1, 'm' : 2}, d_GS = 350E-9, bvK = True)


start_time = time.time()
tau_list2 = SPlt.convergence_tau_plot(twist, TS.Gamma_rot, 150, T = 300, save = True)
print("--- %s seconds ---" % (time.time() - start_time))