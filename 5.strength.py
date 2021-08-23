#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 17:01:45 2021

@author: jiching
"""
import flac
import math
import os,sys
import numpy as np
import pandas as pd
import matplotlib
from math import sqrt
from scipy.special import erf
from matplotlib import cm
import function_savedata as fs
import function_for_flac as f2
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
print('wwwww')
model='w0611'
yy=70
path = '/home/jiching/geoflac/'+model+'/'
path = '/Volumes/My Book/model/'+model+'/'
os.chdir(path)
fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
frame=end

def read_depth(zmesh,x_index,z_index):
    depth=zmesh[x_index,0]-zmesh[x_index,z_index]
    return depth
def nodes_to_elements(xmesh,zmesh,frame):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
x,z=fl.read_mesh(frame)
ele_x,ele_z=nodes_to_elements(x,z,frame)
temp=fl.read_temperature(frame) #on node
phase =fl.read_phase(frame) #on element
visc=fl.read_visc(frame)
strainrate=fl.read_srII(frame)
edot= 1e-15
density = fl.read_density(frame)
g=9.81
pressure=fl.read_pres(frame)
pres = density * g 

v = 10**visc[yy,:]
s = 10**strainrate[yy,:]
p = pressure[yy,:]
dddepth=ele_z[yy,:]
depth=dddepth-dddepth[0]
den=density[yy,:]
visco_strength=2.0 * v * edot
# visco_strength=2.0 * v * s

frico_strength = -p * np.tan(np.pi*(30.0/180.0))
frico_strength2 = -depth  * den * g * np.tan(np.pi*(30.0/180.0))

#-------------------
visco_strength=visco_strength/1e8
frico_strength2=frico_strength2/1e3
#-------------------

fig, ax = plt.subplots(1,1,figsize=(6,10))
applied_strength = np.amin((visco_strength,frico_strength2),axis=0)
ax.plot(visco_strength,depth,'--r',alpha=0.5)
ax.plot(frico_strength2,depth,'--b',alpha=0.5)
ax.plot(applied_strength,depth,'k',lw=3)
ax.set_ylim(-80,0)
ax.set_xlim(0,1000)
ax.set_title('Rock Strength',fontsize=26)
ax.set_xlabel('Strength (MPa)',fontsize=22)
ax.set_ylabel('Depth (km)',fontsize=22)
ax.grid()
# fig.savefig('/home/jiching/geoflac/figure/'+str(model)+'_'+str(yy)+'strength'+'.png')


def half_space_cooling_T(z, Tsurf, Tmantle,  age_in_myrs):
    diffusivity = 1e-6
    myrs2sec = 86400 * 365.2425e6

    T = Tsurf + (Tmantle - Tsurf) * erf(z /
            sqrt(4 * diffusivity * age_in_myrs * myrs2sec) )
    return T


def get_visc(edot, T, n, A, E):
    '''edot: second invariant of strain rate
    T: temperature in Celsius
    n, A, E: viscosity parameters

    return viscosity in Pascal.s
    '''
    R = 8.31448  # gas constant
    pow = 1.0/n - 1
    pow1 = -1.0/n
    visc = 0.25 * (edot**pow) * (0.75*A)**pow1 * np.exp(E / (n * R * (T + 273))) * 1e6
    return visc


def visc_profile(z, T, edot, layerz, nAEs):
    '''Viscosity profile of multi-layers
    z: numpy array of depth (in meters)
    T: array of temperature (in Celsius)
    edot: strain rate (in 1/second)
    layerz: (0, z1, z2, ...) the depth interface of the layers
    nAEs: ( ..., (n, A, E), ...) visc parameters of each layers
    '''

    if layerz[0] != 0:
        print("Error: layerz[0] is not 0", layerz)
    nlayers = len(layerz)
    layerz = tuple(layerz) + (z[-1],)  # deepest depth

    viscp = np.zeros_like(z)
    for i in range(nlayers):
        n, A, E = nAEs[i][:]
        vs = get_visc(edot, T, n, A, E)

        # find depth range of each layer
        z0, z1 = layerz[i], layerz[i+1]
        n0 = (z >= z0).argmax()
        n1 = (z >= z1).argmax()
        #print(i, z0, fz1, n0, n1)

        viscp[n0:n1] = vs[n0:n1]
    return viscp

# if __name__ == "__main__":

layerz = (0, 1e3, 6e3, 16e3)   # 1st elem must be 0
nAEs = ( (3.00, 5.00e+2, 2.00e+5),
         (3.05, 1.25e-1, 3.76e+5),
         (3.00, 7.00e-5, 5.20e+5),
         (3.00, 7.00e+4, 5.20e+5))
# edot = 1e-14  # high strain rate
# edot = 1e-15  # low strain rate

# upper bound of z
deepz = layerz[-1] * 2

z = np.linspace(0, deepz, num=101)
T = half_space_cooling_T(z, 10, 1330, 30)

visc = visc_profile(z, T, edot, layerz, nAEs)
fig2, ax = plt.subplots(1,1,figsize=(6,10))
visco_strength=2.0 * visc * edot
visco_strength=visco_strength
# applied_strength = np.amin((visco_strength,frico_strength2),axis=0)
ax.plot(visco_strength,-z/1000,'--r',alpha=0.5)
ax.plot(frico_strength2,depth,'--b',alpha=0.5)
# ax.plot(applied_strength,depth,'k',lw=3)
ax.set_ylim(-60,0)                                     
ax.set_xlim(0,100000000000)
ax.set_title('Rock Strength',fontsize=26)
ax.set_xlabel('Strength (MPa)',fontsize=22)
ax.set_ylabel('Depth (km)',fontsize=22)
ax.grid()