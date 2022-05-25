#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:26:14 2022

@author: ji-chingchen
"""


import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import math
import flac
import pandas as pd
import gravity as fg
import function_savedata as fs
import function_for_flac as fd
#import flac_interpolate as fi
fig1=1  ## topo, bourger anomaly, free-air anomaly
fig2=1  ## topo, free-air anomaly, model phase plot
bwith = 3

model = sys.argv[1]
frame = int(sys.argv[2])
path = '/scratch2/jiching/03model/'
path = '/home/jiching/geoflac/'
savepath = '/home/jiching/geoflac/data/'
#savepath = '/Users/ji-chingchen/Desktop/data/'
#savepath='D:/model/data/'
figpath = '/home/jiching/geoflac/figure/'

os.chdir(path+model)
fl = flac.Flac()
time=fl.time

def get_gravity(frame):
    px, topo, topomod, fa_gravity, gb_gravity=fg.compute_gravity2(frame)
    px *= 10**-3
    topo *= 10**-3
    topomod *=10**-3
    fa_gravity *= 10**5
    gb_gravity *= 10**5
    return px, topo, topomod, fa_gravity, gb_gravity
px, topo, topomod, fa_gravity, gb_gravity = get_gravity(frame)
fs.save_5txt(str(model)+'topo-grav_'+str(frame), savepath, px, topo, topomod, fa_gravity, gb_gravity)

if fig1: 
    fig, (ax)= plt.subplots(3,1,figsize=(15,12))   
    name=str(model)+'topo-grav_'+str(frame)+'.txt'
    px, topo, topomod, fa_gravity, gb_gravity=np.loadtxt(savepath+name).T
    trench = px[np.argmin(topo)]
    pxx = px[px>trench]
    topo=topo[px>trench]
    fa_gravity=fa_gravity[px>trench]
    gb_gravity=gb_gravity[px>trench]
    px=pxx
    ax[0].plot(px-trench,topo,c="#000080",lw=5,label='model')
    ax[1].plot(px-trench,gb_gravity,c="#000080",lw=5,label='model')
    ax[2].plot(px-trench,fa_gravity,c="#000080",lw=5,label='model')

    # xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
    # with_plot = (xmean>0)*(ztop<-5)
    # xmean = xmean[with_plot]
    # ztop = ztop[with_plot]
    # ax5.plot(xmean,ztop,c='k',lw=3)
        
    data = np.loadtxt(savepath+'Mexico_free-air.txt')
    x,y,qq,fa = data.T
    sx = x[0]
    sy = y[0]
    faxx=np.zeros(len(x))
    for uu in range(1,len(x)):
        faxx[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax[2].plot(faxx,fa,color='green',lw=5,label = 'observation')
    
    data=np.loadtxt(savepath+'Mexico_topo.txt')
    y,topo = data.T
    sx = x[0]
    sy = y[0]
    new_cord=np.zeros(len(x))
    for uu in range(1,len(x)):
        new_cord[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax[0].plot(new_cord,topo/1000,color='green',lw=5,label='observation')
    G=6.6726e-11
    gb=fa-2*np.pi*1800*G*topo*1e5
    ax[1].plot(new_cord,gb,color='green',lw=5,label = 'observation')
    ax[0].legend(fontsize=24)
#================================figure setting================================
    for uu in range(len(ax)):
        ax[uu].tick_params(axis='x', labelsize=25)
        ax[uu].tick_params(axis='y', labelsize=25)
        ax[uu].set_xlim(0,450)
        ax[uu].grid()
        ax[uu].spines['bottom'].set_linewidth(bwith)
        ax[uu].spines['top'].set_linewidth(bwith)
        ax[uu].spines['right'].set_linewidth(bwith)
        ax[uu].spines['left'].set_linewidth(bwith)
    
    ax[0].set_title(model,fontsize=26)
    ax[0].set_ylabel('Topo (km)',fontsize=25)
    ax[1].set_ylabel('Bourger Gravity anomaly (mgal)',fontsize=25)
    ax[2].set_ylabel('Free-air Gravity anomaly (mgal)',fontsize=25)
    ax[-1].set_xlabel('Distance (km)',fontsize=20)

    fig.savefig(figpath+'gravity_vs_topo'+str(model)+'_'+str(frame)+'.png')
    #fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/master03/Seminar/my present/ch0810.png')
#====================================figure2===================================
if fig2:
    fig2, (ax)= plt.subplots(3,1,figsize=(15,12))
    name=str(model)+'topo-grav_'+str(frame)+'.txt'
    px, topo, topomod, fa_gravity, gb_gravity=np.loadtxt(savepath+name).T
    trench = px[np.argmin(topo)]
    pxx = px[px>trench]
    topo=topo[px>trench]
    fa_gravity=fa_gravity[px>trench]
    px=pxx
    ax[0].plot(px-trench,topo,c="#000080",lw=5,label='model')
    ax[1].plot(px-trench,fa_gravity,c="#000080",lw=5,label='model')
    def nodes_to_elements(xmesh,zmesh):
        ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
        ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
        return ele_x, ele_z
    def plot_snapshot(frame):
        x,z = fl.read_mesh(frame)
        xtop,ztop = fd.get_topo(x,z)
        phase = fl.read_phase(frame)
        ele_x,ele_z = nodes_to_elements(x, z)
        temp = fl.read_temperature(frame)
        return x, z, ele_x, ele_z, phase, temp, ztop
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    ax[2].scatter(ele_x-trench,ele_z,c = phase,cmap = phase15,vmax=20,vmin=1)
    cx=ax[2].contour(x-trench,z,temp,cmap = 'rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=1)
    ax[2].clabel(cx, inline=True, fontsize=10,colors='white',fmt="%1.0f")
#================================figure setting================================
    for uu in range(len(ax)):
        ax[uu].tick_params(axis='x', labelsize=25)
        ax[uu].tick_params(axis='y', labelsize=25)
        ax[uu].set_xlim(0,600)
        ax[uu].grid()
        ax[uu].spines['bottom'].set_linewidth(bwith)
        ax[uu].spines['top'].set_linewidth(bwith)
        ax[uu].spines['right'].set_linewidth(bwith)
        ax[uu].spines['left'].set_linewidth(bwith)

    ax[2].set_ylim(-300,30)
    ax[2].set_ylabel('Depth (km)',fontsize=20)
    ax[-1].set_xlabel('Distance (km)',fontsize=20)

    data = np.loadtxt(savepath+'Mexico_free-air.txt')
    x,y,qq,slab = data.T
    sx = x[0]
    sy = y[0]
    new_cord=np.zeros(len(x))
    for uu in range(1,len(x)):
        new_cord[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax[1].plot(new_cord,slab,color='green',lw=5,label = 'observation')
    data=np.loadtxt(savepath+'Mexico_topo.txt')
    y,slab = data.T
    sx = x[0]
    sy = y[0]
    new_cord=np.zeros(len(x))
    for uu in range(1,len(x)):
        new_cord[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax[0].plot(new_cord,slab/1000,color='green',lw=5,label='observation')
    ax[0].legend(fontsize=24)
    fig2.savefig(figpath+'gravity_topo_phase'+str(model)+'_'+str(frame)+'.png')
