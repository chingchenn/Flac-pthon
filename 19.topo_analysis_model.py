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
model_list=['h0401','h0402','h0403','h0404']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
for kk,model in enumerate(model_list):
    xmean,ztop=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_stack_topography.txt').T
    ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=3)

ax2.set_xlim(-width,width)
ax2.set_title("Topography comparation")
ax2.set_ylabel("Bathymetry (m)")
ax2.set_xlabel("Distance relative to trench (km)",fontsize=20)
ax2.legend(fontsize=16)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+'multi_topo_analysis.jpg')
