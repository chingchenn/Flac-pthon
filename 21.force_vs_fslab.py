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
import function_for_flac as fd

#model = sys.argv[1]
model_list=['Chi01','chih0601']#,'chih0602','chih0603']
#model_list=['chimex0601','chimex0602','chimex0603','chimex0604','chimex0605','chimex0606','chimex0607','chimex0608','chimex0609','chimex0610','chimex0611']
model_list=['h0805','h0810','h0808']
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F'] 
savepath='/home/jiching/geoflac/data/'
fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
fig3, (ax4)= plt.subplots(1,1,figsize=(6,8))   
for kk,model in enumerate(model_list):
    path = '/home/jiching/geoflac/data/'+model+'_forc.txt'
    #path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
    #path = '/Volumes/My Book/model/'+model+'/'
    temp1=np.loadtxt(path)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax.scatter(time,forc_l,label=model,color=newcolors[kk],s=2)
    ax2.scatter(time,forc_r,label=model,color=newcolors[kk],s=2)
    movvl = fd.moving_window_smooth(vl,400)
    ax4.scatter(time,movvl*31545741325,label=model,color=newcolors[kk],s=2)
    name=model+'_flatslab_time_len.txt'
#    if os.fstat((savepath+name).fileno()).st_size:
#        data = np.loadtxt(savepath+name, unpack=True)
#    qqaa = np.loadtxt(savepath+name)
#    if len(qqaa) > 0:
#        time,length,depth=np.loadtxt(savepath+name).T
#        ax.axvspan(time[0],time[-1],facecolor=newcolors[kk], alpha=0.2)
    
ax.axhline(y=1.134747173859944351e13,xmin=0,xmax=24,color='r',linestyle='dashed')
ax.set_xlim(0,time[-1])
ax.set_title('oceanic side force',fontsize=16)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid()
ax2.set_xlim(0,time[-1])
ax2.set_title('continental side force',fontsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
ax2.legend()
ax4.set_xlim(0,time[-1])
ax4.set_title('oceanic side velocity',fontsize=16)
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax4.grid()
ax4.set_xlabel('Time (Myr)',fontsize=16)
ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
ax4.legend()
fig.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_forc_flat.png')
fig3.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_vel_TS.png')
