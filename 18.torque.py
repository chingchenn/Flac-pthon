#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 23:04:20 2022

@author: ji-chingchen
"""
import math
import flac
import os,sys
import numpy as np
from scipy import interpolate
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_for_flac as f2
import function_savedata as fs
model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model+'/'
path = '/scratch2/jiching/03model/'+model+'/'
#path = '/Volumes/My Book/model/'+model+'/'
#path = 'F:/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
###================================= parameter setting =======================================
depth1=0
depth2=-200
phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
ecden = 3480-3300
ocden=2800-3300
g=10
thickness=6000
viscosity = 1e24                                                                        #kg/m/s^2 *s
U = 5*3.17*10**-10 
###======================================Time Series==========================================
Torque_G = np.zeros(end)
Torque_H = np.zeros(end)
###===================================find lithisohere========================================
for i in range(20,end):
    x, z = fl.read_mesh(i)
    ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
    ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
    phase = fl.read_phase(i)
    trench_ind = np.argmin(z[:,0]) 
    crust_xc = np.zeros(nex);crust_zc = np.zeros(nex)
    crust_xe = np.zeros(nex);crust_ze = np.zeros(nex)
    for j in range(trench_ind,nex):
        ind_oceanicc = (phase[j,:] == phase_oceanic) 
        ind_oceanice = (phase[j,:] == phase_ecolgite)
        if True in ind_oceanicc:
            crust_xc[j] = np.average(ele_x[j,ind_oceanicc])
            crust_zc[j] = np.average(ele_z[j,ind_oceanicc])
        if True in ind_oceanice:
            crust_xe[j] = np.average(ele_x[j,ind_oceanice])
            crust_ze[j] = np.average(ele_z[j,ind_oceanice])

    ind_within = (crust_zc >= int(depth2)) * (crust_zc < int(depth1))        
    ind_withine = (crust_ze >= int(depth2)) * (crust_ze < int(depth1))

    crust_xmin = np.amin(crust_xc[ind_within]);crust_xmax = np.amax(crust_xc[ind_within])
    crust_zmin = np.amin(crust_zc[ind_within]);crust_zmax = np.amax(crust_zc[ind_within])
    dx = crust_xmax - crust_xmin;dz = crust_zmax - crust_zmin
    anglec= math.degrees(math.atan(dz/dx))*np.pi/180

    crust_xmine = np.amin(crust_xe[ind_withine]);crust_xmaxe = np.amax(crust_xe[ind_withine])
    crust_zmine = np.amin(crust_ze[ind_withine]);crust_zmaxe = np.amax(crust_ze[ind_withine])
    edx = crust_xmaxe - crust_xmine;edz = crust_zmaxe - crust_zmine
    anglee= math.degrees(math.atan(edz/edx))*np.pi/180

    #rc=abs((max(crust_xc)-ele_x[trench_ind,0])/np.cos(anglec))*1e3
    #re=abs((max(crust_xe)-max(crust_xc))/np.cos(anglee))*1e3
    #Torque_G = rc**2*ocden*g*np.sin((90-anglec))-re*ecden*g*np.sin((90-anglee)) # kg*m/s^2
    basalt_length = (max(crust_xc)-ele_x[trench_ind,0])/np.cos(anglec) * 1e3    # m
    eclogite_length = (max(crust_xe)-max(crust_xc))/np.cos(anglee) *1e3          # m
    bc = g*thickness*basalt_length*ocden
    be = g*thickness*eclogite_length*ecden
    Torque_G[i] = 0.5*bc*(basalt_length)**2*np.cos(anglec)+0.5*be*(eclogite_length)**2*np.cos(anglee) # kg*m/s^2

    ## Hydrodynamic                                                                      #m/s
    PA = np.sin(anglec)/((np.pi-anglec)+np.sin(anglec))
    PB = np.sin(anglec)**2/((anglec)**2-np.sin(anglec)**2)
    Torque_H[i]  = 2*viscosity*U*(basalt_length+eclogite_length)*(PA+PB)                       #kg*m/s^2

###============================
fig, (ax)= plt.subplots(1,1,figsize=(13,8))
ax.scatter(fl.time,abs(Torque_G),c='#6A5ACD',label='gravity torque')
ax.scatter(fl.time,Torque_H,c="#D2691E",label='hydrodynamic troque')
ax.legend(fontsize=16)
ax.set_xlabel('Time (Myr)',fontsize=16)
ax.set_ylabel('Force (kg*m/s^2)',fontsize=16)
ax.set_xlim(0, fl.time[-1])
ax.set_title('Torque of '+model,fontsize=20,fontname="Times New Roman")
fig.savefig('/home/jiching/geoflac/figure/'+model+'_torque.png')
fs.save_2txt(model+'torque','/home/jiching/geoflac/data/',Torque_G,Torque_H)
