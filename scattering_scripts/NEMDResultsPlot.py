#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 17:18:41 2019

@author: ramyagurunathan

TBC versus temperature plots
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Qing Hao measured data

'''
Log-log plots comparing to Qing Hao's data
'''
qh_data = np.loadtxt('../datafiles/Si_tbr.csv', delimiter = ',')
qh_data[:,1] = qh_data[:,1]*1e-9
#Data
plt.figure()
plt.scatter(qh_data[:,0], qh_data[:,1], label = 'QH data')
plt.xscale('log')
plt.yscale('log')

plt.loglog(qh_data[:,0], qh_data[0,1]*(qh_data[:,0]/qh_data[0,0])**(-1), label = r'T$^{-1}$')
plt.loglog(qh_data[:,0], qh_data[0,1]*(qh_data[:,0]/qh_data[0,0])**(-1.53), label = r'T$^{-1.53}$')
plt.legend()
plt.xlabel('T (K)')
plt.ylabel(r'$R_K$ $(m^2K/W)$')
plt.ylim([1e-9,1e-6])
plt.savefig('qhTdepedence.pdf', bbox_inches = 'tight')




'''
NEMD Plots for heterointerfaces
'''

NEMD_data = np.loadtxt('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/qinghao_sige_rk/NEMD_het_interface.csv', delimiter = ',')

plt.figure()
NEMD_data[:,1] = 1/NEMD_data[:,1]
plt.scatter(NEMD_data[:,0], NEMD_data[:,1])
plt.xscale('log')
plt.yscale('log')

plt.loglog(NEMD_data[:,0], NEMD_data[0,1]*(NEMD_data[:,0]/NEMD_data[0,0])**(-0.5), label = r'T$^{-0.5}$')
plt.loglog(NEMD_data[:,0], NEMD_data[0,1]*(NEMD_data[:,0]/NEMD_data[0,0])**(-1), label = r'T$^{-1}$')
plt.legend()
plt.xlabel('T (K)')
plt.ylabel(r'$R_K$ $(m^2K/W)$')


'''
Minnich Al on Si Data
'''

AlSi_data = np.loadtxt('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/qinghao_sige_rk/MinnichDataAlSi.csv', delimiter = ",")

plt.figure()
AlSi_data[:,1] = 1/AlSi_data[:,1]
plt.scatter(AlSi_data[:,0], AlSi_data[:,1])
plt.xscale('log')
plt.yscale('log')

plt.loglog(AlSi_data[:,0], AlSi_data[0,1]*(AlSi_data[:,0]/AlSi_data[0,0])**(-0.5), label = r'T$^{-0.5}$')
plt.loglog(AlSi_data[:,0], AlSi_data[0,1]*(AlSi_data[:,0]/AlSi_data[0,0])**(-3), label = r'T$^{-3}$')
plt.legend()
plt.xlabel('T (K)')
plt.ylabel(r'$R_K$ $(m^2K/W)$')


'''
Si-Si Twist Interface Conparison
'''

SiSi_twist = np.loadtxt('/Users/ramyagurunathan/Documents/PhDProjects/The-Grid-Interface/qinghao_sige_rk/Si_Si_twist_datasets.csv', delimiter = ",", skiprows = 2)

angles = ['86.5', '70.4', '15.4', '6.9' ,'3.4']

data = {}

i = 0
for c in range(0,SiSi_twist.shape[1],2):
    data[angles[i]] = [SiSi_twist[:,c], SiSi_twist[:,c+1]]
    i = i+1

plt.figure()
for k,v in data.items():
    plt.scatter(v[0], v[1], label = k)
    plt.xscale('log')
    plt.yscale('log')
    
'''
Fit temperature dependencies
'''
#plt.loglog(data['86.5'][0], data['86.5'][1][0]*(data['86.5'][0]/data['86.5'][0][0])**(-.8), label = r'T$^{-.8}$')
#plt.loglog(data['6.9'][0], data['6.9'][1][0]*(data['6.9'][0]/data['6.9'][0][0])**(-.5), label = r'T$^{-.5}$')

plt.xlabel('T (K)')
plt.ylabel(r'$R_K$ $(m^2K/W)$')   


'''

'''