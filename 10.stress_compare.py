#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:14 2021

@author: ji-chingchen
"""


import os,sys
import numpy as np
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

#model = sys.argv[1]
model_list=['h0410','h0411','h0412','h0413','h0414']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
fig2, (ax3)= plt.subplots(1,1,figsize=(10,8))   
fig3, (ax4)= plt.subplots(1,1,figsize=(10,8))   
for kk,model in enumerate(model_list):
    path = '/home/jiching/geoflac/'+model+'/forc.0'
    #path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
    #path = '/Volumes/My Book/model/'+model+'/'
    temp1=np.loadtxt(path)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax.scatter(time,forc_l,label=model,c=newcolors[kk])
    print(model)
    ax2.scatter(time,forc_r,label=model,c=newcolors[kk])
    ax3.scatter(time,ringforce,label=model,c=newcolors[kk],s=4)
    ax4.plot(time,vl*31545741325,label=model,c=newcolors[kk],lw=2)

ax.set_xlim(0,time[-1])
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid()
ax2.set_xlim(0,time[-1])
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
ax2.legend()
ax3.set_xlim(0,time[-1])
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.grid()
ax3.legend()
ax4.set_xlim(0,time[-1])
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax4.grid()
ax4.set_xlabel('Time (Myr)',fontsize=16)
ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
ax4.legend()
# ax2.plot([time[0],time[-1]],[0,0],'r--')
fig.savefig('/home/jiching/geoflac/figure/forc.png')
fig2.savefig('/home/jiching/geoflac/figure/ringforc.png')
fig3.savefig('/home/jiching/geoflac/figure/vel.png')
