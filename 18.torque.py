#!/usr/bin/env python3
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
# model = str(sys.argv[1])
model='w1260'
path = '/home/jiching/geoflac/'+model+'/'
path = '/Volumes/My Book/model/'+model+'/'
i=97
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
depth1=0
depth2=-120

phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
ecden = 3480-3300
ocden=2800-3300
g=9.8
thickness=6
angle = np.zeros(end)
x, z = fl.read_mesh(i)
ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
phase = fl.read_phase(i)
trench_ind = np.argmin(z[:,0]) 
crust_xc = np.zeros(nex)
crust_zc = np.zeros(nex)
crust_xe = np.zeros(nex)
crust_ze = np.zeros(nex)
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

crust_xmin = np.amin(crust_xc[ind_within])
crust_xmax = np.amax(crust_xc[ind_within])
crust_zmin = np.amin(crust_zc[ind_within])
crust_zmax = np.amax(crust_zc[ind_within])
dx = crust_xmax - crust_xmin
dz = crust_zmax - crust_zmin
anglec= math.degrees(math.atan(dz/dx))
anglee= math.degrees(math.atan((np.amax(crust_ze[ind_withine])-np.amin(crust_ze[ind_withine]))
                               /(np.amax(crust_xe[ind_withine])-np.amin(crust_xe[ind_withine]))))
rc=abs((max(crust_xc)-ele_x[trench_ind,0])/np.cos(anglec*np.pi/180))*1e3
re=abs((max(crust_xe)-max(crust_xc))/np.cos(anglee*np.pi/180))*1e3
torque = rc*ocden*g*np.sin(90-anglec)+re*ecden*g*np.sin(90-anglee) # m^2 ?
fig, (ax)= plt.subplots(1,1,figsize=(13,8))
ax.scatter(crust_xc[ind_within],crust_zc[ind_within])
ax.scatter(crust_xe[ind_withine],crust_ze[ind_withine])
ax.set_aspect('equal')