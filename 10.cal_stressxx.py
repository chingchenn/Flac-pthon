#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:14 2021

@author: ji-chingchen
"""


import flac
import os
import numpy as np
import matplotlib.pyplot as plt

model='b0505m'
path = '/home/jiching/geoflac/data/'
figpath = '/home/jiching/geoflac/figure/'
print('-----plotting flatslab-----')
bwith = 3
name=model+'_flatslab_time_len.txt'
time,length,depth=np.loadtxt(path+name).T
fig2, (ax) = plt.subplots(5,1,figsize=(20,32))
ax[0].scatter(time,depth,c="#000080",s=12)
ax[4].set_xlabel('Time (Myr)',fontsize=16)
ax[0].set_ylabel('depth (km) ',fontsize=16)
ax[1].set_ylabel('Left force (N/m)',fontsize=16)
ax[2].set_ylabel('Right force (N/m) ',fontsize=16)
ax[3].set_ylabel('gravity torque (N*m)',fontsize=16)
ax[4].set_ylabel('suction torque (N*m)',fontsize=16)

for qq in range(len(ax)):
    ax[qq].set_xlim(0,time[-1])
    ax[qq].tick_params(axis='x', labelsize=16)
    ax[qq].tick_params(axis='y', labelsize=16)
    ax[qq].grid()
    ax[qq].spines['bottom'].set_linewidth(bwith)
    ax[qq].spines['top'].set_linewidth(bwith)
    ax[qq].spines['right'].set_linewidth(bwith)
    ax[qq].spines['left'].set_linewidth(bwith)

print('=========== DONE =============')
temp1=np.loadtxt(path+model+"_forc.txt")
nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
ax[1].scatter(time,forc_l,c="#6A5ACD",s=3)
ax[2].scatter(time,forc_r,c = '#6A5ACD',s=3)
temp2=np.loadtxt(path+model+"_torque.txt")
time,Torque_G,Torque_H = temp2.T
ax[3].scatter(time,Torque_G,c="#660000",s=12)
ax[4].scatter(time,Torque_H,c="#660000",s=12)
fig2.savefig(figpath+model+'_slab&force.png')
