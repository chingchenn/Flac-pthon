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
import function_for_flac as fd
import matplotlib.pyplot as plt

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
    xtop,ztop = fd.get_topo(x,z)
    arc_ind,trench_ind=fd.find_trench_index(z)
    trench_x=xtop[trench_ind]
    phase = fl.read_phase(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    temp = fl.read_temperature(frame)
    tren_x = (ele_z[:,0]).argmin()
    ele_x = ele_x-ele_x[tren_x,0]
    return x-trench_x,z,ele_x,ele_z,phase,temp,ztop
#------------------------------------------------------------------------------
x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
fig, (ax)= plt.subplots(1,1,figsize=(12,8))
colors = ["#93CCB1","#8BFF8B","#7158FF","#FF966F","#9F0042",
       "#660000","#524B52","#D14309","#5AB245","#004B00",
       "#008B00","#455E45","#B89FCE","#C97BEA","#525252",
       "#FF0000","#00FF00","#FFFF00","#7158FF"]
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#F87431","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#C97BEA","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
#in_plot=(x>-200)*(x<500)
ax.scatter(ele_x,ele_z,c = phase,cmap = phase15,vmax=19,vmin=1,s=20)
cx=ax.contour(x,z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=1)
#ax.clabel(cx, inline=True, fontsize=10,colors='k',fmt="%1.0f",inline_spacing=0.5)
ax.set_xlim(-200,500)
ax.set_ylim(-300,max(ztop)+2)
ax.set_title(str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
ax.set_ylabel('Depth (km)',fontsize=20)
ax.set_xlabel('Distance from trench (km)',fontsize=20)
ax.set_aspect('equal')
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16 )
fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot.png')
