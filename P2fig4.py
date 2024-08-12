#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 13:08:25 2024

@author: chingchen
"""
import math
import flac
import os,sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib as mpl
from matplotlib import cm
from netCDF4 import Dataset
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation/'


phase_uppercrust = 2
phase_oceanic = 3
phase_mantle1 = 4
phase_schist = 5
phase_mantle2 = 8
phase_serpentinite = 9
phase_sediment = 10
phase_sediment_1 = 11
phase_eclogite = 13
phase_lowercrust = 14
phase_hydratedmantle = 16
phase_oceanic_1 = 17
phase_eclogite_1 = 18

xmin,xmax=250,1000
zmin,zmax=-200,10
model = 'Nazca_aa06'
shift = 320
frame = 51
frame = 11
frame = 190
frame = 126
frame = 200
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
end = 228
bwith=2
fontsize=20
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
              ]

cmap = plt.cm.get_cmap('gist_earth')
# cmap = plt.cm.get_cmap('summer')
fig4, (ax)= plt.subplots(1,1,figsize=(8,12))
fig5, (ax5)= plt.subplots(1,1,figsize=(10,5))
ax.set_title(model,fontsize = 20)
tip_x = np.zeros(end)
tip_z = np.zeros(end)
tttime = np.zeros(end)
x_change_ec = np.zeros(end)
z_change_ec = np.zeros(end)
x_change_serp = np.zeros(end)
z_change_serp = np.zeros(end)
x_change_cho = np.zeros(end)
z_change_cho = np.zeros(end)
x_arc=np.zeros(end)


# xxx_trench = np.max(x_trench)
x, z = fl.read_mesh(150)
chamber_limit = 1e-3
xtt = x[:,0]
ztt = z[:,0]
trench_xx=xtt[np.argmin(ztt)]
# trench = np.zeros(end)
for frame in range(1,end):
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x,z)
    magma = fl.read_fmelt(frame)
    phase = fl.read_phase(frame)
    xt = x[:,0]
    zt = z[:,0]
    # trench_x=xt[np.argmin(zt)]
    # trench[frame]=trench_x
    t = np.zeros(xt.shape)
    t[:] = frame*0.2
    melt_thod = 1e-3
    sxx = fl.read_sxx(frame)*100
    ### topography
    from sklearn.linear_model import LinearRegression
    X = xt
    X = np.reshape(X, (len(X), 1))
    y = zt
    model = LinearRegression()
    model.fit(X, y)
    trend = model.predict(X)
    cb_plot=ax.scatter(x[:,0],t,c=z[:,0]+2,cmap=cmap,vmin=-3, vmax=3,s=20)
    cb_plot2=ax5.scatter(t,x[:,0]-trench_xx,c=2*(z[:,0]-trend),cmap=cmap,vmin=-6, vmax=6,s=10)
    
    # magma 
    if len(ele_x[magma>melt_thod])!=0:
        kk = np.ones(len(ele_x[magma>melt_thod]))*frame*0.2
        ax.scatter(ele_x[magma>melt_thod],kk,c='#FF6347', s=5)
        ax5.scatter(kk,ele_x[magma>melt_thod]-trench_xx,c='#FF6347', s=5)
    # ### arc
    # if True in (phase[:,0]==phase_lowercrust):
    #     kk = np.ones(len(ele_x[:,0][[phase[:,0]==phase_lowercrust]]))*frame*0.2
    #     # ax.scatter(ele_x[:,0][[phase[:,0]==phase_lowercrust]],kk,c='#FF6347', s=5)
    #     # print(ele_x[:,0][[phase[:,0]==phase_lowercrust]])
    if len(sxx[sxx>50])!=0:
        tip_sxx = sxx[sxx>50][-1]
        tip_z[frame] = ele_z[np.where(sxx==tip_sxx)]
        tip_x[frame] = ele_x[np.where(sxx==tip_sxx)]
    tttime[frame] = fl.time[frame]
    
    ### eclogite
    for xx in range(len(ele_x)):
        if True in (phase[xx,:]==phase_eclogite):
            for zz in range(len(ele_x[0])):
                if (phase[xx,zz]==phase_eclogite):
                    x_change_ec[frame]=ele_x[xx,zz]
                    z_change_ec[frame]=-ele_z[xx,zz]
                    break
            break
    
    ### serpentinite
    for xx in range(len(ele_x)-1,0,-1):
        if True in (phase[xx,:]==phase_serpentinite):
            for zz in range(len(ele_x[0])-1,0,-1):
                if (phase[xx,zz]==phase_serpentinite):
                    x_change_serp[frame]=ele_x[xx,zz]
                    z_change_serp[frame]=-ele_z[xx,zz]
                    break
            break
    
    ### chlorite
    for xx in range(len(ele_x)-1,0,-1):
        if True in (phase[xx,:]==phase_hydratedmantle):
            for zz in range(len(ele_x[0])-1,0,-1):
                if (phase[xx,zz]==phase_hydratedmantle):
                    x_change_cho[frame]=ele_x[xx,zz]
                    z_change_cho[frame]=-ele_z[xx,zz]
                    break
            break
ax.scatter(x_change_ec,tttime,c=newcolors[1],s=10,label = 'basalt-eclogite')
ax.scatter(x_change_serp,tttime,c=newcolors[0],s=10,label = 'serpentinite-peridotite')
ax.scatter(x_change_cho,tttime,c=newcolors[3],s=10,label = 'chlorite-peridotite')
ax5.scatter(np.max(ele_x[magma>1.5e-5]),frame*0.2,c='#FF6347',label = 'melt',s = 5)
ax5.scatter(tttime,x_change_ec-trench_xx,c=newcolors[1],s=10)
ax5.scatter(tttime,x_change_serp-trench_xx,c=newcolors[0],s=10)
ax5.scatter(tttime,x_change_cho-trench_xx,c=newcolors[3],s=10)

# ---------------------- plot setting --------------------------
ax.set_xlim(300,xmax)
for aa in [ax,ax5]:
    aa.tick_params(labelsize=fontsize)
    aa.set_ylim(10,40)
    aa.legend(fontsize=fontsize-6)
    aa.set_xlabel('time (Myr)',fontsize=fontsize)
    aa.set_ylabel('distance (km)',fontsize=fontsize)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)

ax_cbin = fig4.add_axes([0.93, 0.125, 0.02, 0.757])
cb = fig4.colorbar(cb_plot,cax=ax_cbin,orientation='vertical')
ax_cbin = fig5.add_axes([0.93, 0.125, 0.02, 0.757])
cb2 = fig5.colorbar(cb_plot2,cax=ax_cbin,orientation='vertical')

for cc in [cb,cb2]:
    cc.set_label(label='bathymetry (km)', size=25)
    cc.ax.yaxis.set_label_position('right')
    cc.ax.tick_params(labelsize=fontsize)



ax5.set_xlim(300,600)
ax5.set_ylim(28,40)
ax5.set_ylim(0,600)
ax5.set_xlim(10,40)
ax5.grid()
# fig5.savefig(figpath+'figure2b_v4.pdf')
# ax5.scatter(trench,tttime,c='black',s=10)