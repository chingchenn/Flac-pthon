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
import function_for_flac as fd
import function_savedata as fs
model = str(sys.argv[1])
#model = 'ch1404'
path = '/home/jiching/geoflac/'+model+'/'
#path = '/scratch2/jiching/03model/'+model+'/'
#path = '/Volumes/SSD500/model/'+model+'/'
# path = '/Volumes/My Book/model/'+model+'/'
#path = 'F:/model/'+model+'/'
savepath = '/scratch2/jiching/data/'
figpath = '/scratch2/jiching/figure/'

###================================= parameter setting =======================================
depth1=0
depth2=-750
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
bet=1
###======================================Time Series==========================================
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
Torque_G = np.zeros(end)
Torque_H = np.zeros(end)
PBB = np.zeros(end)
PAA = np.zeros(end)
anc = np.zeros(end)
ane = np.zeros(end)
ddz = np.zeros(end)
ddx = np.zeros(end)
###======================================Functions==========================================
def trench(end):
    trench_x=np.zeros(end)
    trench_z=np.zeros(end)
    trench_index=np.zeros(end)
    arc_x=np.zeros(end)
    arc_z=np.zeros(end)
    arc_index=np.zeros(end)
    for i in range(1,end):
        x,z = fl.read_mesh(i)
        sx = x[:,0];sz = z[:,0]
        arc_ind,trench_ind=fd.find_trench_index(z)
        trench_index[i]=trench_ind
        arc_index[i]=arc_ind
        trench_x[i]=sx[trench_ind]
        trench_z[i]=sz[trench_ind]
        arc_x[i]=sx[arc_ind]
        arc_z[i]=sz[arc_ind]
    return trench_index,trench_x,trench_z,arc_index,arc_x,arc_z
trench_index,trench_x,trench_z,arc_index,arc_x,arc_z = trench(end)
def find_slab(i):
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)
    xc_ocean = mx[(phase==phase_oceanic)]
    zc_ocean = mz[(phase==phase_oceanic)]
    xe_ocean = mx[(phase==phase_ecolgite)]
    ze_ocean = mz[(phase==phase_ecolgite)]
    start = math.floor(trench_x[i]-50)
    #final = math.floor(np.max(xe_ocean))
    final = math.floor(trench_x[i]+500)
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
        ozc[yy] = np.min(zc_ocean[(xc_ocean>=px1)*(xc_ocean<=xx)])
        oxc[yy] = np.average(xc_ocean[(xc_ocean>=px1)*(xc_ocean<=xx)])
        px1 = xx
    for yy,xx in enumerate(x_grid):
        if len(ze_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])==0:
            continue
        oze[yy] = np.min(ze_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])
        oxe[yy] = np.average(xe_ocean[(xe_ocean>=px2)*(xe_ocean<=xx)])
        px2 = xx
    oxxc=oxc[oxc>start]
    ozc=ozc[oxc>start]
    oxc=oxxc
    oxxe=oxe[oxe>start]
    oze=oze[oxe>start]
    oxe=oxxe
    return oxc,ozc,oxe,oze

###===================================find lithosphere========================================
for i in range(10,end):
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
    if len(zc_ocean[(xc_ocean>=start) *(xc_ocean<=start+bet)])==0:
            continue
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
    anc[i] = anglec
    ddz[i] = dz
    ddx[i] = dx
    basalt_length = (max(oxc)-trench_x[i])/np.cos(anglec) * 1e3
    bc = g*thickness*basalt_length*ocden
    if len(oxe)==0:
        Torque_G[i] = 0.5*bc*(basalt_length)**2*np.cos(anglec)
        continue 
    edx = max(oxe)-min(oxe);edz = max(oze)-min(oze)
    anglee= math.degrees(math.atan(edz/edx))*np.pi/180
    anc[i] = anglec
    ane[i] = anglee
    eclogite_length = (max(oxe)-max(oxc))/np.cos(anglee) *1e3
    be = g*thickness*eclogite_length*ecden
    Torque_G[i] = 0.5*bc*(basalt_length)**2*np.cos(anglec)+0.5*be*(eclogite_length)**2*np.cos(anglee) # kg*m/s^2

    ## Hydrodynamic                                                                  
    PA = np.sin(anglec)/((np.pi-anglec)+np.sin(anglec))
    PB = np.sin(anglec)**2/((anglec)**2-np.sin(anglec)**2)
    PBB[i] = PB
    PAA[i] = PA
    Torque_H[i]  = 2*viscosity*U*(basalt_length+eclogite_length)*(PA+PB)                  	      #kg*m/s^2
fig, (ax)= plt.subplots(1,1,figsize=(10,6))
ax.scatter(fl.time[abs(Torque_G)>0],Torque_G[abs(Torque_G)>0],c='#6A5ACD',label='gravity torque')
ax.scatter(fl.time[Torque_H>0],Torque_H[Torque_H>0],c="#D2691E",label='hydrodynamic troque')
ax.legend(fontsize=16,loc='upper right')
ax.set_xlabel('Time (Myr)',fontsize=16)
ax.set_ylabel('Force (kg*m/s^2)',fontsize=16)
ax.set_xlim(0, fl.time[-1])
# ax.set_yscale('log')
ax.set_title('Torque of '+model,fontsize=20,fontname="Times New Roman")
fig.savefig(figpath+model+'_torque.png')
fs.save_3txt(model+'_torque',savepath,fl.time,Torque_G,Torque_H)


###==================== find crust by another way, (not a good way) =====================
# for i in range(1,end):
#     x, z = fl.read_mesh(i)
#     ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
#     ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
#     phase = fl.read_phase(i)
#     trench_ind = np.argmin(z[:,0]) 
#     crust_xc = np.zeros(nex);crust_zc = np.zeros(nex)
#     crust_xe = np.zeros(nex);crust_ze = np.zeros(nex)
#     for j in range(trench_ind,nex):
#         ind_oceanicc = (phase[j,:] == phase_oceanic) 
#         ind_oceanice = (phase[j,:] == phase_ecolgite)
#         if True in ind_oceanicc:
#             crust_xc[j] = np.average(ele_x[j,ind_oceanicc])
#             crust_zc[j] = np.average(ele_z[j,ind_oceanicc])
#         if True in ind_oceanice:
#             crust_xe[j] = np.average(ele_x[j,ind_oceanice])
#             crust_ze[j] = np.average(ele_z[j,ind_oceanice])

#     ind_within = (crust_zc >= int(depth2)) * (crust_zc < int(depth1))        
#     ind_withine = (crust_ze >= int(depth2)) * (crust_ze < int(depth1))

#     crust_xmin = np.amin(crust_xc[ind_within]);crust_xmax = np.amax(crust_xc[ind_within])
#     crust_zmin = np.amin(crust_zc[ind_within]);crust_zmax = np.amax(crust_zc[ind_within])
#     dx = crust_xmax - crust_xmin;dz = crust_zmax - crust_zmin
#     anglec= math.degrees(math.atan(dz/dx))*np.pi/180

#     if len(crust_xe[ind_withine])==0:
#         continue
#     crust_xmine = np.amin(crust_xe[ind_withine]);crust_xmaxe = np.amax(crust_xe[ind_withine])
#     crust_zmine = np.amin(crust_ze[ind_withine]);crust_zmaxe = np.amax(crust_ze[ind_withine])
#     edx = crust_xmaxe - crust_xmine;edz = crust_zmaxe - crust_zmine
#     anglee= math.degrees(math.atan(edz/edx))*np.pi/180
#     anc[i] = anglec
#     ane[i] = anglee
#     ddz[i] = dz
#     ddx[i] = dx
#     #Gravity
#     #rc=abs((max(crust_xc)-ele_x[trench_ind,0])/np.cos(anglec))*1e3
#     #re=abs((max(crust_xe)-max(crust_xc))/np.cos(anglee))*1e3
#     #Torque_G = rc**2*ocden*g*np.sin((90-anglec))-re*ecden*g*np.sin((90-anglee)) # kg*m/s^2
#     basalt_length = (max(crust_xc)-ele_x[trench_ind,0])/np.cos(anglec) * 1e3
#     eclogite_length = (max(crust_xe)-max(crust_xc))/np.cos(anglee) *1e3
#     bc = g*thickness*basalt_length*ocden
#     be = g*thickness*eclogite_length*ecden
#     Torque_G[i] = 0.5*bc*(basalt_length)**2*np.cos(anglec)+0.5*be*(eclogite_length)**2*np.cos(anglee) # kg*m/s^2

#     ## Hydrodynamic                                                                  
#     PA = np.sin(anglec)/((np.pi-anglec)+np.sin(anglec))
#     PB = np.sin(anglec)**2/((anglec)**2-np.sin(anglec)**2)
#     PBB[i] = PB
#     PAA[i] = PA
#     Torque_H[i]  = 2*viscosity*U*(basalt_length+eclogite_length)*(PA+PB)                  	      #kg*m/s^2

# ###============================
# fig, (ax)= plt.subplots(1,1,figsize=(10,6))
# ax.scatter(fl.time[abs(Torque_G)>0],Torque_G[abs(Torque_G)>0],c='#6A5ACD',label='gravity torque')
# ax.scatter(fl.time[Torque_H>0],Torque_H[Torque_H>0],c="#D2691E",label='hydrodynamic troque')
# ax.legend(fontsize=16,loc='upper right')
# ax.set_xlabel('Time (Myr)',fontsize=16)
# ax.set_ylabel('Force (kg*m/s^2)',fontsize=16)
# ax.set_xlim(0, fl.time[-1])
# # ax.set_yscale('log')
# ax.set_title('Torque of '+model,fontsize=20,fontname="Times New Roman")
# # fig.savefig('/home/jiching/geoflac/figure/'+model+'_torque.png')
# # fs.save_2txt(model+'torque','/home/jiching/geoflac/data/',Torque_G,Torque_H)
