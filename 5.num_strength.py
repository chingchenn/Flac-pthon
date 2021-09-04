#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 23:30:45 2021

@author: ji-ching chen
"""
import flac
import os
import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt

model='w1011'
yy=70
path = '/home/jiching/geoflac/'+model+'/'
os.chdir(path)
fl = flac.Flac()
end = fl.nrec
frame=2
layerz = (0, 18e3, 30e3)   # 1st elem must be 0
Dfc = ((2800,30,4e7),
       (2900,30,4e7),
       (3300,30,4e7))
nAEs = ( (3.05, 1.25e-1, 2.76e+5),
       (3.05, 1.25e-1, 3.76e+5),
       (3.00, 7.00e+4, 5.20e+5))
edot = 1e-14  # high strain rate
edot = 1e-15  # low strain rate
deepz = layerz[-1] * 3
z = np.linspace(0, deepz, num=1000)

#------------------------------------------------------------------------------
#numerical solution of plastic stress
def nodes_to_elements(xmesh,zmesh,frame):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
x,z=fl.read_mesh(frame)
ele_x,ele_z=nodes_to_elements(x,z,frame)
density = fl.read_density(frame)
g=9.81
pressure=fl.read_pres(frame)
p = pressure[yy,:]
dddepth=ele_z[yy,:]
depth=dddepth-dddepth[0]
den=density[yy,:]
frico_strength= -depth * 1000 * den * g * np.tan(np.pi*(30.0/180.0))+4e7 
#------------------------------------------------------------------------------
# equation soluiton of plastic stress and viscosity
layerz = (0, 18e3, 30e3)   # 1st elem must be 0
Dfc = ((2800,30,4e7),
       (2900,30,4e7),
       (3300,30,4e7))
nAEs = ( (3.05, 1.25e-1, 2.76e+5),
       (3.05, 1.25e-1, 3.76e+5),
       (3.00, 7.00e+4, 5.20e+5))
edot = 1e-14  # high strain rate
edot = 1e-15  # low strain rate
deepz = layerz[-1] * 3
z = np.linspace(0, deepz, num=1000)
con_T = f2.continental_geothermal_T(z,20, 6,45)
visc = f2.visc_profile(z, con_T, edot, layerz, nAEs)
visco_strength=visc* edot *2 #Pa
#------------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(6,10))
applied_strength = np.amin((visco_strength,frico_strength),axis=0)
ax.plot(visco_strength/1e6,-z/1000,'--r',alpha=0.5)
ax.plot(frico_strength/1e6,-z/1000,'--b',alpha=0.5)
ax.plot(applied_strength/1e6,-z/1000,'k',lw=3)
ax.set_ylim(-deepz/1000 +1,0)                                     
ax.set_xlim(0,1e3)
ax.set_title('Rock Strength',fontsize=26)
ax.set_xlabel('Strength (MPa)',fontsize=22)
ax.set_ylabel('Depth (km)',fontsize=22)
ax.grid()
fig.savefig('/home/jiching/geoflac/figure/'+'strength_profile_numerical'+'.png')