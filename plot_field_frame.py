#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 10:53:35 2023

@author: chingchen
"""

import flac
import os
import numpy as np
import matplotlib
from netCDF4 import Dataset
import matplotlib.pyplot as plt

### make phase.grd file in 34.phase_interpolate.py

plt.rcParams["font.family"] = "Arial"
#---------------------------------- DO WHAT -----------------------------------
# Model
Cocos           = 1
Nazca           = 0
### pdf or png
png             = 0
pdf             = 0

### plot
shot_12         = 0
shot_5          = 0  ### phase, viscosity, density, dynamic pressure, srII
plot_phase      = 1
plot_viscosity  = 0
plot_density    = 0
plot_pressure   = 0
plot_sxx        = 0
plot_sII        = 0
plot_srII       = 0
plot_phase_topo = 0
shot_3          = 0 ### phase, viscosity, dynamic pressure


#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/Users/chingchen/Desktop/GSA2023/mp4file/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Nazca/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
# figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/AGU/POSTER/Poster_figure/'

if Cocos:
    xmin,xmax=500,900
    zmin,zmax=-120,10
    model = 'Ref_Cocos'
    shift = 550
    # model = 'Cocos11'
    # xmin,xmax=400,900
    # zmin,zmax=-150,10
# elif Nazca:
#     xmin,xmax=250,1000
#     zmin,zmax=-200,10
#     model = 'Ref_Nazca'
elif Nazca:
    xmin,xmax=250,1000
    zmin,zmax=-200,10
    model = 'Nazca_aa06'
    shift = 320
    # model = 'Ref_Nazca'
frame = 51
frame = 11
frame = 190
frame = 130
# frame = 151
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
colors2=[
 '#C98F49', '#92C0DF', '#2553C7', '#FFFFFF', '#6495ED',
 '#2E8B57', '#524B52', '#9A32CD', '#6B8E23','#D4DBF4',
 '#D8BFD8','#999999','#F2C85B','#92C0DF','#999999',
 '#4CC552','#999999','#999999','#999999','#999999']
phase8= matplotlib.colors.ListedColormap(colors2)
end = 201

for frame in range(2,end,2):
    fig, (ax)= plt.subplots(1,1,figsize=(17,7))
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    #file=model+'_frame'+str(frame)+'_phase.grd'
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=2)
    #--------------------- melting plotting -------------------------
    melt=fl.read_fmelt(frame)
    ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='#CD5C5C')
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(labelsize=23)
    ax.set_ylabel('depth (km)',fontsize=26)
    ax.set_xlabel('distance (km)',fontsize=26)
    if Nazca:
        xmajor_ticks = np.linspace(250,1000,num=7)
        ymajor_ticks = np.linspace(200,0,num=5)
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
        ymajor_ticks = np.linspace(150,0,num=4)
    ax.set_xticks(xmajor_ticks)
    ax.set_yticks(ymajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')