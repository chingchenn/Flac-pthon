#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 4 17:26:14 2022

@author: ji-chingchen
"""


import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

model = sys.argv[1]
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F'] 
savepath='/home/jiching/geoflac/data/'
fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(12,8))   
path = '/home/jiching/geoflac/'+model+'/forc.0'
temp1=np.loadtxt(path)
nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
ax.scatter(time,forc_l,label=model,color='#2F4F4F',s=5)
name=model+'_flatslab_time_len.txt'
time,length,depth=np.loadtxt(savepath+name).T
ax.axvspan(time[0],time[-1],facecolor='#52254F', alpha=0.2)
ax3.plot(time,length,c="#000080",lw=3)
time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
ax2.scatter(time[melt>0],xmelt[melt>0],c=melt[melt>0],cmap='gray',s=45)
    
ax.axhline(y=1.134747173859944351e13,xmin=0,xmax=24,color='r',linestyle='dashed')
#ax.set_xlim(0,time[-1])
ax.set_title('oceanic side force',fontsize=16)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid()
#ax2.set_xlim(0,time[-1])
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
ax.set_xlim(0,30)
ax2.set_xlim(0,30)
ax3.set_xlim(0,30)
ax2.set_ylim(200,600)
ax2.set_ylabel('Distance (km)',fontsize=20)
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.set_xlabel('Time (Myr)',fontsize=16)
ax3.grid()
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)
ax3.spines['bottom'].set_linewidth(bwith)
ax3.spines['top'].set_linewidth(bwith)
ax3.spines['right'].set_linewidth(bwith)
ax3.spines['left'].set_linewidth(bwith)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_forc_flat_melt.png')
