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

width=800
depth1=-5
depth2=-140
fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
model_list=['h0420','h0421','h0422','h0423','h0424','h0425','h0426']
model_list=['h0427','h0428','h0429','h0430','h0431','h0432','h0433']
model_list=['h0409','h0408','h0405','h0406','h0407']
model_list=['chih0605','chih0614','chih0615','chih0616','chih0617']
newcolors = ['#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
fig, (ax)= plt.subplots(1,1,figsize=(10,7))
for kk,model in enumerate(model_list):
    xmean,ztop=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_stack_slab.txt').T
    ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=3)
    df = pd.read_csv('/home/jiching/geoflac/data/plate_dip_of_'+model+'.csv')
    ax.plot(df.time[df.angle>0],df.angle[df.angle>0],c=newcolors[kk],lw=2,label =model)
    print(newcolors[kk])


ax2.set_xlim(-10,500)
ax2.set_title("slab comparation")
ax2.set_ylabel("Depth (km)")
ax2.set_xlabel("Distance relative to trench (km)",fontsize=16)
ax2.legend(fontsize=16)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.png')
ax.set_title('Angle Variation of '+str(model),fontsize=24)
ax.set_xlabel('Time (Myr)',fontsize=20)
ax.set_ylabel('Angel ($^\circ$) from '+str(-depth1)+' to '+str(-depth2)+' depth',fontsize=20)
ax.grid()
ax.legend(fontsize=16)
fig.savefig('/home/jiching/geoflac/'+'figure/'+model_list[0]+'_'+model_list[-1]+'_dip.jpg')
