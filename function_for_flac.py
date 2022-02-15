#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 13:28:16 2021

@author: jiching
"""


import sys
import math
import pandas as pd
import numpy as np
from math import sqrt
from scipy.special import erf

sys.path.append('/home/jiching/geoflac/util')

def get_topo(xmesh,zmesh):
    xtop=xmesh[:,0]
    ztop=zmesh[:,0]
    return xtop,ztop

def find_trench_index(z):
    zz = z[:,0]
    if zz.all()==0:
        imax,i=0,0
    else:
        imax = zz.argmax()
        i = zz[:imax].argmin()
    return imax,i

def read_depth(z_array,x_index,z_index):
    depth=z_array[x_index,0]-z_array[x_index,z_index]
    return depth

def read_area(xmesh,zmesh,x_index,z_index):
    x1 = xmesh[x_index,z_index]
    y1 = zmesh[x_index,z_index]
    x2 = xmesh[x_index,z_index+1]
    y2 = zmesh[x_index,z_index+1]
    x3 = xmesh[x_index+1,z_index+1]
    y3 = zmesh[x_index+1,z_index+1]
    x4 = xmesh[x_index+1,z_index]
    y4 = zmesh[x_index+1,z_index]
    area1 = ((x1-x2)*(y3-y2))-((x3-x2)*(y1-y2))
    area2 = ((x1-x4)*(y3-y4))-((x3-x4)*(y1-y4))    
    area = (abs(area1)+abs(area2))*0.5           
    return area
"""
def read_time(start_vts,model_steps):
    timestep=[0]
    for step in range(start_vts,model_steps+1):
        timestep.append(fl.time[step])
    # timestep=np.array(timestep)
    return timestep
"""

#find melt element
def melt_element(xmesh,zmesh,frame,mm):
    melt_xele=[]
    melt_zele=[]
    melt_number=[]
    for xx in range(len(mm)):
        for zz in range(len(mm[0])):
            if mm[xx,zz] != 0:
                melt_xele.append(xx)
                melt_zele.append(zz)
                melt_number.append(mm[xx,zz])
    return melt_xele,melt_zele,melt_number

def chamber_element(xmesh,zmesh,frame,mm):
    chamber_xele=[]
    chamber_zele=[]
    chamber_number=[]
    for xx in range(len(mm)):
        for zz in range(len(mm[0])):
            if mm[xx,zz] != 0:
                chamber_xele.append(xx)
                chamber_zele.append(zz)
                chamber_number.append(mm[xx,zz])
    return chamber_xele,chamber_zele,chamber_number

def moving_window_smooth(A,window_width):
    MM = np.zeros(len(A))    
    for kk in range(0,len(A)):
        if kk>(len(A)-window_width):
            MM[kk] = A[kk]
        else:
            MM[kk] = sum(A[kk:kk+window_width])/window_width
    return MM
def half_space_cooling_T(z, Tsurf, Tmantle,  age_in_myrs):
    diffusivity = 1e-6
    myrs2sec = 86400 * 365.2425e6

    T = Tsurf + (Tmantle - Tsurf) * erf(z /
            sqrt(4 * diffusivity * age_in_myrs * myrs2sec) )
    return T

def continental_geothermal_T3(z,cond1,cond2,depth):
    T = np.zeros_like(z)
    for kk,zz in enumerate(z):
        if zz/1000 > depth:
            T[kk]= depth * cond1 + (zz/1000-depth) * cond2 
        else:T[kk] = cond1 * zz/1000  
    T[T>1330]=1330
    return T
def continental_geothermal_T4(z, Tsurf, Tmantle, depth):
    T = np.zeros_like(z)
    T = Tsurf + (Tmantle - Tsurf)/depth * z/1000
    T[T>1330]=1330
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

def plastic_stress(z,layerz,Dfc,g=9.81):

    nlayers = len(layerz)
    layerz = tuple(layerz) + (z[-1],)  

    plast = np.zeros_like(z)
    for i in range(nlayers):
        den, fric1, coh1 = Dfc[i][:]
        ps= z * den * g * np.tan(np.pi*(fric1/180.0))+coh1 
        z0, z1 = layerz[i], layerz[i+1]
        n0 = (z >= z0).argmax()
        n1 = (z >= z1).argmax()
        plast[n0:n1] = ps[n0:n1]
    return plast

def getDistance(latA, lonA, latB, lonB):
    ra = 6378140  
    rb = 6356755  
    flatten = (ra - rb) / ra  # Partial rate of the earth
    # change angle to radians
    radLatA = math.radians(latA)
    radLonA = math.radians(lonA)
    radLatB = math.radians(latB)
    radLonB = math.radians(lonB)

    pA = math.atan(rb / ra * math.tan(radLatA))
    pB = math.atan(rb / ra * math.tan(radLatB))
    x = math.acos(math.sin(pA) * math.sin(pB) + math.cos(pA) * math.cos(pB) * math.cos(radLonA - radLonB))
    c1 = (math.sin(x) - x) * (math.sin(pA) + math.sin(pB)) ** 2 / math.cos(x / 2) ** 2
    c2 = (math.sin(x) + x) * (math.sin(pA) - math.sin(pB)) ** 2 / math.sin(x / 2) ** 2
    dr = flatten / 8 * (c1 - c2)
    distance = ra * (x + dr)
    distance = round(distance / 1000, 4)
    return distance
def phase_pro(phasenumber):
    ff=pd.read_csv('phase.csv',index_col='phase')
    density,alfa,bata,n,A,E,rl,rm,plas1,plas2,fric1,fric2,coh1,coh2,dilat1,dilat2,cond,cp,Ts,Tl,Tk,fk=ff.iloc[phasenumber-1]
    return density,alfa,bata,n,A,E,rl,rm,plas1,plas2,fric1,fric2,coh1,coh2,dilat1,dilat2,cond,cp,Ts,Tl,Tk,fk