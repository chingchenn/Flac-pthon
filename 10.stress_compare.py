#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:14 2021

@author: ji-chingchen
"""


import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

#model = sys.argv[1]
#model_list=['chimex0601','chimex0602','chimex0603']
model_list=['chimex0601','chimex0602','chimex0603','chimex0604','chimex0605','chimex0606','chimex0607','chimex0608','chimex0609','chimex0610','chimex0611']
#model_list=['h0401','h0402','h0403','h0404']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = ['#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F'] 
fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
fig2, (ax3)= plt.subplots(1,1,figsize=(10,8))   
fig3, (ax4)= plt.subplots(1,1,figsize=(10,8))   
for kk,model in enumerate(model_list):
#    if kk==0:
    path = '/home/jiching/geoflac/'+model+'/forc.0'
#    else:
#    path = '/scratch2/jiching/03model/'+model+'/forc.0'
    #path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
    #path = '/Volumes/My Book/model/'+model+'/'
    temp1=np.loadtxt(path)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax.scatter(time,forc_l,label=model,color=newcolors[kk],s=2)
    ax2.scatter(time,forc_r,label=model,color=newcolors[kk],s=2)
    ax3.scatter(time,ringforce,label=model,color=newcolors[kk],s=2)
    ax4.scatter(time,vl*31545741325,label=model,color=newcolors[kk],s=2)

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
ax3.set_xlim(0,time[-1])
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.grid()
ax3.legend()
ax4.set_xlim(0,time[-1])
ax4.set_title('oceanic side velocity',fontsize=16)
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax4.grid()
ax4.set_xlabel('Time (Myr)',fontsize=16)
ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
ax4.legend()
# ax2.plot([time[0],time[-1]],[0,0],'r--')
fig.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_forc.png')
fig2.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_ringforc.png')
fig3.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_vel.png')
