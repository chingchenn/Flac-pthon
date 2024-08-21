#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 22:32:34 2023

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


model_list=['Nazca_aa06','Nazca_ab05','Nazca_ab03','Nazca_ab01']#,'Nazca_ab08']
label_list=['reference','11 km ridge','3580 kg/m3','16 km ridge','3630 kg/m3']
color_list=['#CD5C5C','#FF8C00','#228B22','#4682B4','#8B008B','#','#']
bwith=3

fig, (ax1)= plt.subplots(1,1,figsize=(15,8))
#### time 
i = 175    ## 10 Myr 
frame1 = i*0.2
for kk,model in enumerate(model_list):
    xx,zz,cx = np.loadtxt(savepath+model+'_'+str(frame1)+'_final_slab.txt').T
    xx=xx[zz<0]
    cx=cx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    #ax1.plot(xx,-zz,color=color_list[kk],lw=5)
    ax1.plot(xx,-zz,color=color_list[kk],lw=5,label=label_list[kk])
ax1.legend(fontsize=20)


# fig1, (ax2)= plt.subplots(1,1,figsize=(15,8))
# for kk,model in enumerate(model_list):
#     xx,zz,mx = np.loadtxt(savepath+model+'_'+str(int(frame1))+'_final_moho_slab.txt').T
#     xx=xx[zz<0]
#     mx=mx[zz<0]
#     zz=zz[zz<0]
#     zz = fd.moving_window_smooth(zz,3)
# #    ax2.plot(xx,-zz,color=color_list[kk],lw=5)
#     ax2.plot(mx,-zz,color=color_list[kk],lw=5)

for ax in [ax1]:#,ax2]:
    ax.grid()
    ax.set_xlim(0,700)
    ax.set_ylim(150,0)
    ax.set_aspect('equal')
    ax.set_ylabel('depth (km)',fontsize=26)
    ax.set_xlabel('distance (km)',fontsize=26)
    ax.tick_params(labelsize=26,width = 2)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
# #fig.savefig(figpath+'fig5av6.pdf')
# #fig1.savefig(figpath+'fig5av5.pdf')


model_list=['Nazca_aa06','Nazca_v2_06','Nazca_v2_ab03','Nazca_v2_ab01','Nazca_v2_ab04']
label_list=['reference ori','reference new','3580 kg/m3','16 km ridge','3380 kg/m3']

fig3, (ax3)= plt.subplots(1,1,figsize=(15,8))

for kk,model in enumerate(model_list):
    xx,zz,cx = np.loadtxt(savepath+model+'_'+str(frame1)+'_final_slab.txt').T
    xx=xx[zz<0]
    cx=cx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    #ax1.plot(xx,-zz,color=color_list[kk],lw=5)
    ax3.plot(xx,-zz,color=color_list[kk],lw=5,label=label_list[kk])
ax3.legend(fontsize=20)
ax3.grid()
ax3.set_xlim(0,700)
ax3.set_ylim(150,0)
ax3.set_aspect('equal')
ax3.set_ylabel('depth (km)',fontsize=26)
ax3.set_xlabel('distance (km)',fontsize=26)
ax3.tick_params(labelsize=26,width = 2)
for axis in ['top','bottom','left','right']:
    ax3.spines[axis].set_linewidth(bwith)