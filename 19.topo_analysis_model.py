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

width=700
fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
fig3, (ax3,ax10,ax11) = plt.subplots(3,1,figsize=(15,8))
fig4, (ax4,ax5,ax6) = plt.subplots(3,1,figsize=(15,8))
fig5,(ax7,ax8,ax9) = plt.subplots(3,1,figsize=(15,8))
#model_list=['chih0609','chih0613','chih0617']
model_list=['chih0605','chih0614','chih0615','chih0616','chih0617']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
for kk,model in enumerate(model_list):
    xmean,ztop=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_stack_topography.txt').T
    ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=3)
    time,distance,xtrench,xarc=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_location_data.txt').T
    time,arc_basin_distance,trench_basin_distance,basin_width=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_dis_data.txt').T
    time,basin_depth,arc_height,trench_height=np.loadtxt('/home/jiching/geoflac/data/'+str(model)+'_depth_data.txt').T
    ax3.plot(time,distance,lw=2,label=model,color=newcolors[kk])
    ax10.plot(time,xtrench,lw=2,label=model,color=newcolors[kk])
    ax11.plot(time,xarc,lw=2,label=model,color=newcolors[kk])
    ax4.plot(time,arc_basin_distance,lw=2,label=model,color=newcolors[kk])
    ax5.plot(time,trench_basin_distance,lw=2,label=model,color=newcolors[kk])
    ax6.plot(time,basin_width,lw=2,label=model,color=newcolors[kk])
    ax7.plot(time,basin_depth,lw=2,label=model,color=newcolors[kk])
    ax8.plot(time,arc_height,lw=2,label=model,color=newcolors[kk])
    ax9.plot(time,trench_height,lw=2,label=model,color=newcolors[kk])


ax2.set_xlim(-width,width)
ax2.set_title("Topography comparation")
ax2.set_ylabel("Bathymetry (km)")
ax2.set_xlabel("Distance relative to trench (km)",fontsize=16)
ax2.legend(fontsize=16)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+'multi_topo_analysis'+model_list[0]+'_'+model_list[-1]+'.jpg')
ax3.grid()
ax3.legend(fontsize=16)
ax3.set_ylim(0,time[-1])
ax3.set_title("ATdis comparation")
ax3.set_ylabel("Time (Ma)",fontsize=16)
ax3.set_xlabel("Distance differences (km)",fontsize=16)
ax10.grid()
ax11.grid()
ax3.set_xlim(0,time[-1])
ax3.set_title('arc trench distance',fontsize=12)
ax10.set_xlim(0,time[-1])
ax10.set_title('trench location (km)',fontsize=12)
ax11.set_xlim(0,time[-1])
ax11.set_title('arc location (km)',fontsize=12)
fig3.savefig('/home/jiching/geoflac/figure'+'/'+'multi_distance'+model_list[0]+'_'+model_list[-1]+'.jpg')
ax4.set_xlim(0,time[-1])
ax4.grid()
ax4.set_title('arc basin distance',fontsize=12)
ax5.set_xlim(0,time[-1])
ax5.grid()
ax5.set_title('trench basin distance',fontsize=12)
ax6.set_xlim(0,time[-1])
ax6.grid()
ax6.set_title('basin width',fontsize=12)
ax7.set_xlim(0,time[-1])
ax7.grid()
ax7.set_title(model+'\n basin depth',fontsize=12)
ax8.set_xlim(0,time[-1])
ax8.grid()
ax8.set_title('arc height',fontsize=12)
ax9.set_xlim(0,time[-1])
ax9.grid()
ax9.set_title('trench height',fontsize=12)
fig4.savefig('/home/jiching/geoflac/figure'+'/'+'multi_distance'+model_list[0]+'_'+model_list[-1]+'_xdirection.png')
fig5.savefig('/home/jiching/geoflac/figure'+'/'+'multi_distance'+model_list[0]+'_'+model_list[-1]+'_zdirction.png')
