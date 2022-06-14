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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_for_flac as f2
import function_savedata as fs
model_list = ['ch1404','ch1406']

color=['#2F4F4F','#4682B4','#CD5C5C','#708090',
      '#AE6378','#282130','#7E9680','#24788F',
      '#849DAB','#EA5E51','#35838D','#4198B9',
      '#414F67','#97795D','#6B0D47','#A80359',
      '#52254F','r'] 
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
viscosity = 1e24
U = 5*3.17*10**-10
bet=2
fig, (ax)= plt.subplots(2,1,figsize=(13,8))
for mm,model in enumerate(model_list):
    path = '/home/jiching/geoflac/'+model+'/'
    path = '/scratch2/jiching/03model/'+model+'/'
    #path = '/Volumes/My Book/model/'+model+'/'
    #path = 'F:/model/'+model+'/'
    os.chdir(path)
    fl = flac.Flac();end = fl.nrec
    nex = fl.nx - 1;nez = fl.nz - 1 
    ###======================================Time Series==========================================
    Torque_G = np.zeros(end)
    Torque_H = np.zeros(end)
    ###===================================find lithisohere========================================
    for i in range(2,end):
        x, z = fl.read_mesh(i)
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)
        ## In this code, we considered the marker phase, not the element phase
        trench_ind = np.argmin(z[:,0])
        x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]
        xc_ocean = mx[(phase==phase_oceanic)]
        zc_ocean = mz[(phase==phase_oceanic)]
        xe_ocean = mx[(phase==phase_ecolgite)]
        ze_ocean = mz[(phase==phase_ecolgite)]
        if z_trench> -2 or min(zc_ocean)>-50:
            continue
        start = math.floor(x_trench)
        final = math.floor(np.max(xe_ocean))
        x_grid = np.arange(start,final,bet)
        oxc = np.zeros(len(x_grid))
        ozc = np.zeros(len(x_grid))
        oxe = np.zeros(len(x_grid))
        oze = np.zeros(len(x_grid))
        px1 = start-bet
        px2 = start-bet
        #find initial basalt depth to remove the weage basalt
        kk=np.max(zc_ocean[(xc_ocean>=start) *(xc_ocean<=start+bet)])
        xc_ocean = xc_ocean[zc_ocean<kk]
        zc_ocean = zc_ocean[zc_ocean<kk]
        # interplate to the grid length "bet"
        for yy,xx in enumerate(x_grid):
            if len(zc_ocean[(xc_ocean>=px1)*(xc_ocean<=xx)])==0:
                continue    
            ozc[yy] = np.average(zc_ocean[(xc_ocean>=px1)*(xc_ocean<=xx)])
            oxc[yy] = np.average(xc_ocean[(xc_ocean>=px1)*(xc_ocean<=xx)])
            px1 = xx
        for yy,xx in enumerate(x_grid):
            if len(ze_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])==0:
                continue
            oze[yy] = np.average(ze_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])
            oxe[yy] = np.average(xe_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])
            px2 = xx
        oxxc=oxc[oxc>start]
        ozc=ozc[oxc>start]
        oxc=oxxc
        oxxe=oxe[oxe>start]
        oze=oze[oxe>start]
        oxe=oxxe
        dx = max(oxc)-min(oxc);dz = max(ozc)-min(ozc)
        anglec= math.degrees(math.atan(dz/dx))*np.pi/180
        edx = max(oxe)-min(oxe);edz = max(oze)-min(oze)
        anglee= math.degrees(math.atan(edz/edx))*np.pi/180
        #Gravity
        basalt_length = (max(oxc)-x_trench)/np.cos(anglec) * 1e3
        eclogite_length = (max(oxe)-max(oxc))/np.cos(anglee) *1e3
        bc = g*thickness*basalt_length*ocden
        be = g*thickness*eclogite_length*ecden
        Torque_G[i] = 0.5*bc*(basalt_length)**2*np.cos(anglec)+0.5*be*(eclogite_length)**2*np.cos(anglee) # kg*m/s^2

        ## Hydrodynamic                                                                  
        PA = np.sin(anglec)/((np.pi-anglec)+np.sin(anglec))
        PB = np.sin(anglec)**2/((anglec)**2-np.sin(anglec)**2)
        Torque_H[i]  = 2*viscosity*U*(basalt_length+eclogite_length)*(PA+PB)                  	      #kg*m/s^2

    ###============================
    ax[0].scatter(fl.time[abs(Torque_G)>0],abs(Torque_G)[abs(Torque_G)>0],c=color[mm],label='gravity torque of '+model)
    ax[1].scatter(fl.time[abs(Torque_G)>0],Torque_H[abs(Torque_G)>0],c=color[mm],label='hydrodynamic troque of '+model)
ax[0].legend(fontsize=16)
ax[1].set_xlabel('Time (Myr)',fontsize=16)
ax[0].set_title('Torque of '+model,fontsize=20)
for qq in range(len(ax)):
    ax[qq].set_ylabel('Force (kg*m/s^2)',fontsize=16)
    ax[qq].set_xlim(0, fl.time[-1])
    ax[qq].set_ylim(-0.5e23,3e23)

fig.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_torque.png')
fs.save_2txt(model+'torque','/home/jiching/geoflac/data/',Torque_G,Torque_H)
