#!/usr/bin/env python
# 2022 feb.23
import flac
import sys,os
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs

model = sys.argv[1]
path = '/home/jiching/geoflac/'+model+'/'
#path='/scratch2/jiching/03model/'+model+'/'
#path='/scratch2/jiching/22winter/'+model+'/'
#path = 'F:/model/'+model+'/'
os.chdir(path)
    
fl = flac.Flac();end = fl.nrec
rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))
cmap = plt.cm.get_cmap('gist_earth')
savepath = '/home/jiching/geoflac/data/'
###==========================Time Series===============================
distance = np.zeros(end) #x
xtrench = np.zeros(end) #x
xarc = np.zeros(end) #x
arc_basin_distance = np.zeros(end) #x
trench_basin_distance = np.zeros(end) #x
basin_width = np.zeros(end) #x
basin_depth = np.zeros(end) #z
arc_height = np.zeros(end) #z
trench_height = np.zeros(end) #z
###======================== Stacked Topography ========================
fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
topo1 = 0;xmean = 0;ictime = 20
width = 600
for i in range(1,end):
    x, z = fl.read_mesh(i)
    xt = x[:,0]
    zt = z[:,0]
    t = np.zeros(xt.shape)
    t[:] = i*0.2
    x_trench = xt[np.argmin(zt)]
    xtrench[i]=x_trench
    x_arc = xt[np.argmax(zt[xt<x_trench+200])]
    xarc[i]=x_arc
    z_trench = zt[np.argmin(zt)]
    z_arc = zt[np.argmax(zt)]
    within_plot = (xt>x_trench-width) * (xt<x_trench+width)
    if i >= end-ictime:
        topo1 += zt
        xmean += (xt-x_trench)
    ax2.plot(xt[within_plot]-x_trench,zt[within_plot],c=newcolors[i]) # Stacked Topography
    distance[i] = x_arc-x_trench 
    if distance[i]>600 or distance[i]<0:
        distance[i]=0
    z_mid_basin = zt[np.argmin(zt[xt>x_arc])]
    x_mid_basin = xt[np.argmin(zt[xt>x_arc])]
    arc_basin_distance[i] = x_mid_basin-x_arc
    trench_basin_distance[i] = x_mid_basin-x_trench
    basin_depth[i] = z_mid_basin
    zb = z_mid_basin/2
    xbw = xt[zt<z_mid_basin]
    if len(xbw)==0:
        continue
    basin_width[i]=(xbw[-1]-xbw[0])
    arc_height[i] = z_arc
    trench_height[i] = z_trench
    
ax2.plot((xmean[within_plot]/ictime), (topo1[within_plot]/ictime),c='purple')
ax2.set_xlim(-width,width)
ax2.set_title(str(model)+"Topography")
ax2.set_ylabel("Bathymetry (km)")
ax2.set_xlabel("Distance relative to trench (km)")
fs.save_2txt(str(model)+'_stack_topography',savepath,xmean[within_plot]/ictime,topo1[within_plot]/ictime)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo_analysis.jpg')
###================ distance between trench and arc peak ================
fig3, (ax3,ax10,ax11) = plt.subplots(3,1,figsize=(15,8))
ax3.plot(fl.time,distance,c='teal',lw=2)
ax10.plot(fl.time,xtrench,c='teal',lw=2)
ax11.plot(fl.time,xarc,c='teal',lw=2)
ax3.grid()
ax10.grid()
ax11.grid()
ax3.set_xlim(0,fl.time[-1])
ax3.set_title('arc trench distance',fontsize=12)
ax10.set_xlim(0,fl.time[-1])
ax10.set_title('trench location (km)',fontsize=12)
ax11.set_xlim(0,fl.time[-1])
ax11.set_title('arc location (km)',fontsize=12)
fs.save_4txt(str(model)+'_location_data',savepath,fl.time,distance,xtrench,xarc)
fig3.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_distance.jpg')
###============================== distance  ==============================
fig4, (ax4,ax5,ax6) = plt.subplots(3,1,figsize=(15,8))
ax4.plot(fl.time,arc_basin_distance,c='teal',lw=2)
ax5.plot(fl.time,trench_basin_distance,c='teal',lw=2)
ax6.plot(fl.time,basin_width,c='teal',lw=2)
ax4.set_xlim(0,fl.time[-1])
ax4.grid()
ax4.set_title('arc basin distance',fontsize=12)
ax5.set_xlim(0,fl.time[-1])
ax5.grid()
ax5.set_title('trench basin distance',fontsize=12)
ax6.set_xlim(0,fl.time[-1])
ax6.grid()
ax6.set_title('basin width',fontsize=12)
fs.save_4txt(str(model)+'_dis_data',savepath,fl.time,arc_basin_distance,trench_basin_distance,basin_width)
fig4.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_xdirection.png')
fig5,(ax7,ax8,ax9) = plt.subplots(3,1,figsize=(15,8))
ax7.plot(fl.time,basin_depth,c='teal',lw=2)
ax8.plot(fl.time,arc_height,c='teal',lw=2)
ax9.plot(fl.time,trench_height,c='teal',lw=2)
ax7.set_xlim(0,fl.time[-1])
ax7.grid()
ax7.set_title(model+'\n basin depth',fontsize=12)
ax8.set_xlim(0,fl.time[-1])
ax8.grid()
ax8.set_title('arc height',fontsize=12)
ax9.set_xlim(0,fl.time[-1])
ax9.grid()
ax9.set_title('trench height',fontsize=12)
fs.save_4txt(str(model)+'_depth_data',savepath,fl.time,basin_depth,arc_height,trench_height)
fig5.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_zdirction.png')
