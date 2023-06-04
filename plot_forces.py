#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:18:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import matplotlib
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.figsize"] = (10,12)
model = sys.argv[1]
#model = 'b0601m'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/03model/'
#path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T

bwith = 3
###----------------------- Slab sinking force with time-------------------------------


fs.save_5txt(model+'_forces','/home/jiching/geoflac/data/',fl.time,fsb,ft,fsu,ratio)
fig, (ax)= plt.subplots(1,1,figsize=(10,6))
sb = fd.moving_window_smooth(fsb[fsb>0],8)
tt = fd.moving_window_smooth(fsu[fsu>0],8)
ax.plot(fl.time[fsb>0],sb,c='#c06c84',label='slab pull (N/m)',lw=4)
ax.plot(fl.time[fsu>0],tt,c="#355c7d",label='suction force (N/m)',lw=4)
#ax.scatter(fl.time[fsb>0],fsb[fsb>0],c='#c06c84',label='slab pull (N/m)')
#ax.scatter(fl.time[ft>0],ft[ft>0],c="#355c7d",label='traction force (N/m)')
ax.legend(fontsize=16,loc='upper left')
#================================figure setting================================
ax.set_xlabel('Time (Myr)',fontsize=16)
ax.set_ylabel('Force (N/m)',fontsize=16)
ax.set_xlim(0, fl.time[-1])
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid()
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
#ax.set_yscale('log')
ax.set_title('Forces of '+model,fontsize=20)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force.png')
fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
ratio_f = fd.moving_window_smooth(ratio[ratio>0],5)
ax2.plot(fl.time[ratio>0],ratio_f,c="#355c7d",label='ratio of these forces)',lw=4)
ax2.set_xlabel('Time (Myr)',fontsize=16)
ax2.set_xlim(0, fl.time[-1])
ax2.set_ylim(0, 30)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)
fig2.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force_ratio.png')
