#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 11:04:51 2022

@author: ji-chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as f2
import matplotlib.pyplot as plt

#---------------------------------- DO WHAT -----------------------------------
shot = 1
gravity = 0
pressure = 0
gravity_plot = 0
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
sys.path.append("/home/jiching/geoflac/util")

model = sys.argv[1]
frame = int(sys.argv[2])
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
    xtop,ztop = f2.get_topo(x,z)
    phase = fl.read_phase(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    temp = fl.read_temperature(frame)
    return x,z,ele_x,ele_z,phase,temp,ztop
def get_gravity(frame):
    px, topo, topomod, fa_gravity, gb_gravity=fg.compute_gravity2(frame)
    px *= 10**-3
    topo *= 10**-3
    topomod *=10**-3
    fa_gravity *= 10**5
    gb_gravity *= 10**5
    return px, topo, topomod, fa_gravity, gb_gravity
def get_pressure(frame):
    x,z = fl.read_mesh(frame)
    xtop,ztop = f2.get_topo(x,z)
    ele_x,ele_z = nodes_to_elements(x, z)
    pre=fl.read_pres(frame)
    onepre=-pre.flatten()
    a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(onepre-fit).reshape(len(pre),len(pre[0]))*100
    return ele_x,ele_z,dypre,ztop
#------------------------------------------------------------------------------
## Creat Data
if gravity:
    if not os.path.exists(savepath+model):
        os.makedirs(savepath+model)
    px, topo, topomod, fa_gravity, gb_gravity = get_gravity(frame)
    fs.save_5array('topo-grav.'+str(frame), savepath+model, px, topo, topomod, fa_gravity,
                     gb_gravity, 'disX', 'topo', 'topomod', 'free-air', 'bourger')
#------------------------------------------------------------------------------
if shot:
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    fig, (ax)= plt.subplots(1,1,figsize=(12,8))
    colors = ["#93CCB1","#8BFF8B","#7158FF","#FF966F","#9F0042",
           "#660000","#524B52","#D14309","#5AB245","#004B00",
           "#008B00","#455E45","#B89FCE","#C97BEA","#525252",
           "#FF0000","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    ax.scatter(ele_x,ele_z,c = phase,cmap = phase15,vmax=20,vmin=1)
    cx=ax.contour(x,z,temp,cmap = 'rainbow',levels = np.arange(100,1500,200),linewidths=1)
    ax.clabel(cx, inline=True, fontsize=10,colors='k')
    ax.set_xlim(0,ele_x[-1,-1])
    ax.set_ylim(ele_z[-1,-1],max(ztop)+2)
    ax.set_title(str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax.set_aspect('equal')
    fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot.png')
if pressure:
    fig, (ax)= plt.subplots(1,1,figsize=(12,8))
    ele_x,ele_z,dypre,ztop = get_pressure(frame)
    cc = plt.cm.get_cmap('RdYlBu')
    cb_plot=ax.scatter(ele_x,ele_z,c=dypre,cmap=cc,vmin=-200, vmax=200)
    ax_cbin = fig.add_axes([0.63, 0.35, 0.23, 0.03])
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    tick_font_size = 10
    cb.ax.tick_params(labelsize=tick_font_size)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax_cbin.set_title('MPa',fontsize=20)
    ax.set_aspect('equal')
    ax.set_xlim(0,max(ele_x[:,0]))
    ax.set_ylim(min(ele_z[0,:]),max(ztop)+2)
    ax.set_title('Pressure '+str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
    fig.savefig(figpath+model+'frame_'+str(frame)+'_dynamic_pressure.pdf')
    
