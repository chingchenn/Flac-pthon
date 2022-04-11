#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 13:45:24 2022

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
from numpy import unravel_index

model = str(sys.argv[1])
figpath='/home/jiching/geoflac/figure/'
path = '/home/jiching/geoflac/'+model+'/'
#path = '/Volumes/My Book/model/'+model+'/'
savepath = '/Users/ji-chingchen/Desktop/data/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1

def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
depth = 0
colors = ["#ffffff","#f0a224"]
phase15= matplotlib.colors.ListedColormap(colors)
fig, (ax)= plt.subplots(1,1,figsize=(10,12))
time=[];ph=[];xx=[]
for step in range(end):
    x, z = fl.read_mesh(step+1)
    phase=fl.read_phase(step+1)
    ele_x, ele_z = nodes_to_elements(x,z)
    xt = ele_x[:,0]
    zt = ele_z[:,0]
    pp = np.zeros(xt.shape)
    t = np.zeros(zt.shape)
    t[:]=fl.time[step]
    for gg in range(len(ele_z)):
        if phase[gg,depth]==14:
            pp[gg]=1
        else:
            pp[gg]=0        
    ax.scatter(xt,t,c=pp,cmap=phase15)
    

ax.set_ylabel("Time (Ma)")
ax.set_xlabel("Distance (km)")
ax.set_title(str(model)+" Phase")
# ax.set_ylim(0,t[-1][-1])
# ax.set_xlim(xt[0][0],xt[-1][-1])
# cb_plot1 = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=phase15,vmin=1, vmax=18)
# ax_cbin = fig.add_axes([0.27, 0.03, 0.23, 0.03])
# cb = fig.colorbar(cb_plot1,cax=ax_cbin,orientation='horizontal')
# ax_cbin.set_title('Phase')
fig.savefig(figpath+model+'_phase.png')
print('=========== DONE =============')