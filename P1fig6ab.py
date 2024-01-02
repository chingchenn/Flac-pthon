#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 14:55:13 2023

@author: chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
from numpy import unravel_index
import matplotlib.colors as mcolors

path = '/home/jiching/geoflac/'
path='/Users/chingchen/Desktop/model/'
#savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
savepath='/Users/chingchen/Desktop/data/'
#figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
figpath='/Users/chingchen/Desktop/figure/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'
plt.rcParams["font.family"] = "Helvetica"


model_list=['Nazca_aa06','Nazca_a0910','Nazca_a0706','Nazca_a0707','Nazca_a0708']
label_list=['5 km','2.5 km','7.5 km','10 km','15 km']
color_list=['#CD5C5C','#FF8C00','#228B22','#4682B4','#8B008B','#','#']
model = model_list[0]
os.chdir(path+model)
fl = flac.Flac()
bwith=2
fontsize=26

fig, (ax1)= plt.subplots(1,1,figsize=(15,8))

#### time 
i = 150    ## 10 Myr 
frame1 = i*0.2
for kk,model in enumerate(model_list):
    xx,zz,cx = np.loadtxt(savepath+model+'_'+str(frame1)+'_final_slab.txt').T
    xx=xx[zz<0]
    cx=cx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    #ax1.plot(cx,-zz,color=color_list[kk],lw=5)
    ax1.plot(xx,-zz,color=color_list[kk],lw=4,label=label_list[kk])



fig1, (ax2)= plt.subplots(1,1,figsize=(15,3))
for kk,model in enumerate(model_list):
    ## PLOT model1 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    tt = fd.moving_window_smooth(fsu,5)/1e19
    ax2.plot(time,tt,c=color_list[kk],label=label_list[kk],lw=4)
    
    # ## PLOT model2 Torque 
    # time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    # tt = fd.moving_window_smooth(fsu,1)/1e19
    # ax2.plot(time,tt,c=color_list[kk],label=label2,lw=4)
    
    # ## PLOT model3 Torque 
    # time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    # tt = fd.moving_window_smooth(fsu,1)/1e19
    # ax2.plot(time,tt,c=color_list[kk],label=label3,lw=4)
    
    # ## PLOT model4 Torque 
    # time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    # tt = fd.moving_window_smooth(fsu,1)/1e19
    # ax2.plot(time,tt,c=color_list[kk],label=label4,lw=4)
#================================figure setting================================
ymajor_ticks = [0,3]
ax2.set_xlim(0,40)
ax2.set_ylim(-1,3)
ax2.set_yticks(ymajor_ticks)
ax2.set_ylabel('suction torque (10$^{19}$ N$\cdot$m/m)',fontsize=fontsize-4)
#ax2.legend(fontsize=fontsize-6,loc='upper left')
ax2.set_xlabel('time (Myr)',fontsize=fontsize)
ymajor_ticks = np.linspace(120,0,num=4)
ax1.set_yticks(ymajor_ticks)
ax1.set_xlim(0,600)
ax1.set_ylim(120,0)
ax1.set_aspect('equal')
ax1.set_yticks(ymajor_ticks)
#ax1.legend(fontsize=20)
ax1.set_ylabel('depth (km)',fontsize=fontsize)
ax1.set_xlabel('distance (km)',fontsize=fontsize)
ax2.set_ylabel('suction torque \n (10$^{19}$ N$\cdot$m/m)',fontsize=fontsize)
ax2.set_xlabel('time (Myr)',fontsize=fontsize)
for ax in [ax1,ax2]:
    ax.grid()
    ax.tick_params(labelsize=fontsize,width = 2)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
fig.savefig(figpath+'fig6av1.pdf')
fig1.savefig(figpath+'fig6bv1.pdf')