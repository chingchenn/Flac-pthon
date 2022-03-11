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
#model = 'h0133'
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
distance = np.zeros(end)
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
    x_arc = xt[np.argmax(zt)]
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
    
ax2.plot((xmean[within_plot]/ictime), (topo1[within_plot]/ictime),c='purple')
ax2.set_xlim(-width,width)
ax2.set_title(str(model)+"Topography")
ax2.set_ylabel("Bathymetry (km)")
ax2.set_xlabel("Distance relative to trench (km)")
fs.save_2txt(str(model)+'_stack_topography',savepath,xmean[within_plot]/ictime,topo1[within_plot]/ictime)
fig2.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo_analysis.jpg')
###================ distance between trench and arc peak ================
fig3, (ax3) = plt.subplots(1,1,figsize=(10,6))
ax3.plot(distance,fl.time,c='teal',lw=4)
ax3.grid()
ax3.set_ylim(0,fl.time[-1])
fs.save_2txt(str(model)+'_dis_trench_arc',savepath,fl.time,distance)
fig3.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_distance.jpg')
