#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:15:39 2022

@author: chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
#import gravity as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
#import flac_interpolate as fi

plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- DO WHAT -----------------------------------
# Model
Cocos           = 1
Nazca           = 0
### pdf or png
png             = 0
pdf             = 0

### plot
shot            = 1
shot_zoomin     = 1
pressure        = 0
viscosity       = 0
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/'

# model = sys.argv[1]
# frame = int(sys.argv[2])

if Cocos:
    xmin,xmax=400,800
    zmin,zmax=150,0
    model = 'Ref_Cocos'
if Nazca:
    xmin,xmax=150,600
    zmin,zmax=150,0
    model = 'Nazca_a0702'
frame = 200
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

#------------------------------------------------------------------------------
def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
def plot_snapshot(frame):
    x,z = fl.read_mesh(frame)
    xtop,ztop = fd.get_topo(x,z)
    phase = fl.read_phase(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    temp = fl.read_temperature(frame)
    return x, z, ele_x, ele_z, phase, temp, ztop
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
#------------------------------------------------------------------------------
if shot:
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    fig, (ax)= plt.subplots(1,1,figsize=(12,8))
    ax.pcolormesh(x,-z,phase,cmap = phase15,vmax=20,vmin=1,shading = 'auto')
    cx=ax.contour(x,-z,temp,cmap = 'rainbow',levels =[200,400,600,800,1000,1200],linewidths=1)
    ax.contour(x,-z,temp,levels =[1300],linewidths=2,colors = '#F08080',linestyles='dashed')
    ax.clabel(cx, inline=True, fontsize=10,colors='white',fmt="%1.0f")
    ax.set_xlim(0,1200)
    ax.set_ylim(300,-20)
    # ax.set_title(str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax.set_aspect('equal')
    bwith = 3
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot.png')
    if pdf:
        if Cocos:
            fig.savefig(figpath+'Ref_Cocos/'+model+'frame_'+str(frame)+'_snapshot.pdf')
        if Nazca:
            fig.savefig(figpath+'Ref_Nazca/'+model+'frame_'+str(frame)+'_snapshot.pdf')
if shot_zoomin:
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    fig, (ax)= plt.subplots(1,1,figsize=(12,8))
    ax.pcolormesh(x,-z,phase,cmap = phase15,vmax=20,vmin=1,shading = 'auto')
    cx=ax.contour(x,-z,temp,cmap = 'rainbow',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax.contour(x,-z,temp,levels =[1300],linewidths=2,colors = '#F08080',linestyles='dashed')
    # ax.clabel(cx, inline=True, fontsize=10,colors='white',fmt="%1.0f")
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)
    # ax.set_title(str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax.set_aspect('equal')
    bwith = 3
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_zoomin.png')
    if pdf:
        if Cocos:
            fig.savefig(figpath+'Ref_Cocos/'+model+'frame_'+str(frame)+'_snapshot_zoomin.pdf') 
        if Nazca:
            fig.savefig(figpath+'Ref_Nazca/'+model+'frame_'+str(frame)+'_snapshot_zoomin.pdf') 
