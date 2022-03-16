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
model_list=['h0420','h0421','h0422','h0423','h0424','h0425','h0426']
model_list=['h0427','h0428','h0429','h0430','h0431','h0432','h0433']
model_list=['h0409','h0408','h0405','h0406','h0407']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
for kk,model in enumerate(model_list):
    xmean,ztop=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_stack_slab.txt').T
    ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=3)

ax2.set_xlim(-10,500)
ax2.set_title("slab comparation")
ax2.set_ylabel("Depth (km)")
ax2.set_xlabel("Distance relative to trench (km)",fontsize=16)
ax2.legend(fontsize=16)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.png')
