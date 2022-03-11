#!/usr/bin/env python
import flac
import sys,os
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs

width=400
fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
fig3, (ax3) = plt.subplots(1,1,figsize=(8,6))
model_list=['h0401','h0402','h0403','h0404']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
for kk,model in enumerate(model_list):
    xmean,ztop=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_stack_topography.txt').T
    ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=3)
    time,distance=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_dis_trench_arc.txt').T
    ax3.plot(distance,time,c=newcolors[kk],label=model,lw=3)

ax2.set_xlim(-width,width)
ax2.set_title("Topography comparation")
ax2.set_ylabel("Bathymetry (km)")
ax2.set_xlabel("Distance relative to trench (km)",fontsize=16)
ax2.legend(fontsize=16)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+'multi_topo_analysis.jpg')
ax3.grid()
ax3.legend(fontsize=16)
ax3.set_ylim(0,time[-1])
ax3.set_title("ATdis comparation")
ax3.set_ylabel("Time (Ma)",fontsize=16)
ax3.set_xlabel("Distance differences (km)",fontsize=16)
fig3.savefig('/home/jiching/geoflac/figure'+'/'+'multi_distance.jpg')
