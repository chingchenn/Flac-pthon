#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 15:17:29 2024

@author: chingchen
"""


import flac
import os
import numpy as np
# import pandas as pd
#import gravity as fg
import matplotlib
import matplotlib as mpl
#matplotlib.use('Agg')
from matplotlib import cm
# from netCDF4 import Dataset
# import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
#import flac_interpolate as fi
plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation/'

xmin,xmax = 250,1000
zmin,zmax = -200,10
frame1 = 30
frame2 = 60
frame3 = 120
frame4 = 178

bwith = 3
fontsize=30
labelsize=30

model='Nazca_aa06'
# model='Nazca_v2_06'
os.chdir(path+model)
fl = flac.Flac()
end = 200
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
cmap = 'jet'

time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
rainbow = cm.get_cmap('gray_r',end)
meltcolor = cm.get_cmap('turbo',end)
newcolors = rainbow(np.linspace(0, 1, end))
time_color = meltcolor(np.linspace(0,1,end))

fig8, (ax2,ax5)= plt.subplots(1,2,figsize=(20,8),gridspec_kw={'width_ratios':[2,0.05]})
xxx_trench = np.max(x_trench)
x, z = fl.read_mesh(150)
chamber_limit = 1e-3
xtt = x[:,0]
ztt = z[:,0]
trench_x=xtt[np.argmin(ztt)]
for i in range(1,end):
    x, z = fl.read_mesh(i)
    ele_x, ele_z = flac.elem_coord(x,z)
    magma_chamber = fl.read_fmagma(i) 
    melt = fl.read_fmelt(i) * 100
    ax2.scatter(ele_x[magma_chamber>chamber_limit],-ele_z[magma_chamber>chamber_limit],c=time_color[i],cmap =cmap,s=10,vmin=0,vmax=40)
for ax in [ax2]:
    ax.grid()
    ax.set_xlim(300,1000)
    ax.set_ylim(150,0)
    ax.set_ylabel('depth (km)',fontsize=fontsize)
    ax.tick_params(labelsize=fontsize )
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.set_aspect('equal')
    ymajor_ticks = np.linspace(150,0,num=4)
    ax.set_yticks(ymajor_ticks)
for frame in [10.0, 20.0, 30.0, 40.0]:
    xx,zz,xt = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
    xx=xx[zz<0]
    zz=zz[zz<0]
    ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
    print(int(frame)*5-1)
    
ax2.set_xlabel('distance (km)',fontsize=fontsize)
norm = mpl.colors.Normalize(vmin=0,vmax=40)
cc1 = mpl.colorbar.ColorbarBase(ax5,cmap=cmap,norm = norm, orientation='vertical')
cc1.set_label(label='time (Myr)', size=fontsize)
cc1.ax.tick_params(labelsize=fontsize)

fig8.savefig(figpath+'figure3_color.pdf')

# fig10, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(20,11),gridspec_kw={'height_ratios':[1.5,2,0.15]})
# xxx_trench = np.max(x_trench)
# chamber_limit = 1e-3
# for i in range(1,end):
#     x, z = fl.read_mesh(i)
#     ele_x, ele_z = flac.elem_coord(x,z)
#     magma_chamber = fl.read_fmagma(i) 
#     melt = fl.read_fmelt(i) * 100
#     ax2.scatter(ele_x[magma_chamber>chamber_limit]-trench_x,-ele_z[magma_chamber>chamber_limit],color=time_color[i],cmap =cmap,s=10,vmin=0,vmax=40)
# for ax in [ax2,ax1]:
#     ax.grid()
#     ax.set_xlim(300,1000)
#     ax.tick_params(labelsize=fontsize )
#     for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(bwith)
#     ymajor_ticks = np.linspace(150,0,num=4)
#     ax.set_yticks(ymajor_ticks)
# for frame in [10.0, 20.0, 30.0, 40.0]:
#     xx,zz,xt = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
#     xx=xx[zz<0]
#     zz=zz[zz<0]
#     ax2.plot(xx,-zz,color=time_color[int(frame)*5-1],lw=3)
#     print(int(frame)*5-1)

# ax2.set_xlim(0,700)
# ax1.set_xlim(0,700)
# ax1.set_ylim(-6,3)
# ymajor_ticks = np.linspace(-6,3,num=4)
# ax1.set_yticks(ymajor_ticks)

# from sklearn.linear_model import LinearRegression
# X = xtt
# X = np.reshape(X, (len(X), 1))
# y = ztt
# model = LinearRegression()
# model.fit(X, y)
# trend = model.predict(X)
# # ax.plot(xtt[:150]-trench_x,ztt[:150],c='k')
# ax.plot(xtt-trench_x,(ztt-trend)*2,c='b',lw=3)
# ax2.set_aspect('equal')
# ax2.set_ylim(150,0)
# ax1.set_ylabel('topography (km)',fontsize=fontsize) 
# ax2.set_ylabel('depth (km)',fontsize=fontsize)   
# ax2.set_xlabel('distance (km)',fontsize=fontsize)
# norm = mpl.colors.Normalize(vmin=0,vmax=40)
# cc1 = mpl.colorbar.ColorbarBase(ax5,cmap=cmap,norm = norm, orientation='horizontal')
# cc1.set_label(label='time (Myr)', size=fontsize)
# cc1.ax.tick_params(labelsize=fontsize)

# # for axis in ['top','bottom','left','right']:
# #     ax.spines[axis].set_linewidth(bwith)
# # fig10.savefig(figpath+'figure2a_v5.pdf')