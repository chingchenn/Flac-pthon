#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 4 17:26:14 2022

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
fig1=1
fig2=1

model = sys.argv[1]
frame=150
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#6B0D47','#A80359','#52254F'] 
path='/scratch2/jiching/03model/'
savepath='/home/jiching/geoflac/data/'
#savepath = '/Users/ji-chingchen/Desktop/data/'
figpath = '/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

if fig1: 
    fig, (ax3,ax4,ax5)= plt.subplots(3,1,figsize=(15,12))   
    name=str(model)+'topo-grav_'+str(frame)+'.txt'
    px, topo, topomod, fa_gravity, gb_gravity=np.loadtxt(savepath+name).T
    trench = px[np.argmin(topo)]
    pxx = px[px>trench]
    topo=topo[px>trench]
    fa_gravity=fa_gravity[px>trench]
    px=pxx
    ax3.plot(px-trench,topo,c="#000080",lw=3)
    ax4.plot(px-trench,fa_gravity,c="#2F4F4F",lw=3)
    xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
    with_plot = (xmean>0)*(ztop<-5)
    xmean = xmean[with_plot]
    ztop = ztop[with_plot]
    ax5.plot(xmean,ztop,c='k',lw=3)
    #================================figure setting================================
    ax3.set_title(model,fontsize=26)
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax5.tick_params(axis='x', labelsize=16)
    ax5.tick_params(axis='y', labelsize=16)
    
    ax3.set_xlim(0,600)
    ax4.set_xlim(0,600)
    ax5.set_xlim(0,600)
    
    #ax3.set_xlabel('Distance (km)',fontsize=20)
    ax3.set_ylabel('Topo (km)',fontsize=20)
    ax4.set_ylabel('mgal',fontsize=20)
    #ax4.set_xlabel('Distance (km)',fontsize=20)
    ax5.set_ylabel('Depth (km)',fontsize=20)
    #ax5.set_xlabel('Distance (km)',fontsize=20)
    

    ax3.grid()
    ax4.grid()
    ax5.grid()
    
    bwith = 3
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    ax5.spines['bottom'].set_linewidth(bwith)
    ax5.spines['top'].set_linewidth(bwith)
    ax5.spines['right'].set_linewidth(bwith)
    ax5.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+'gravity_vs_topo'+str(model)+'_'+str(frame)+'.png')
    #fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/master03/Seminar/my present/ch0810.png')
#====================================figure2===================================
if fig2:
    fig2, (ax3,ax4,ax5)= plt.subplots(3,1,figsize=(15,12))
    name=str(model)+'topo-grav_'+str(frame)+'.txt'
    px, topo, topomod, fa_gravity, gb_gravity=np.loadtxt(savepath+name).T
    trench = px[np.argmin(topo)]
    pxx = px[px>trench]
    topo=topo[px>trench]
    fa_gravity=fa_gravity[px>trench]
    px=pxx
    ax3.plot(px-trench,topo,c="#000080",lw=5,label='model')
    ax4.plot(px-trench,fa_gravity,c="#000080",lw=5,label='model')
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
    ax5.scatter(ele_x-trench,ele_z,c = phase,cmap = phase15,vmax=20,vmin=1)
    cx=ax5.contour(x-trench,z,temp,cmap = 'rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=1)
    ax5.clabel(cx, inline=True, fontsize=10,colors='white',fmt="%1.0f")
    ax3.set_xlim(0,600)
    ax4.set_xlim(0,600)
    ax5.set_xlim(0,600)
    ax5.set_ylim(-300,30)
    ax5.set_ylabel('Depth (km)',fontsize=20)
    ax5.set_xlabel('Distance (km)',fontsize=20)
    bwith = 3
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax3.tick_params(axis='x', labelsize=16 )
    ax3.tick_params(axis='y', labelsize=16 )
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    ax4.tick_params(axis='x', labelsize=16 )
    ax4.tick_params(axis='y', labelsize=16 )
    ax5.spines['bottom'].set_linewidth(bwith)
    ax5.spines['top'].set_linewidth(bwith)
    ax5.spines['right'].set_linewidth(bwith)
    ax5.spines['left'].set_linewidth(bwith)
    ax5.tick_params(axis='x', labelsize=16 )
    ax5.tick_params(axis='y', labelsize=16 )
    data = np.loadtxt(savepath+'Mexico_free-air.txt')
    x,y,qq,slab = data.T
    sx = x[0]
    sy = y[0]
    new_cord=np.zeros(len(x))
    for uu in range(1,len(x)):
        new_cord[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax4.plot(new_cord,slab,color='green',lw=5,label = 'observation')
    data=np.loadtxt(savepath+'Mexico_topo.txt')
    y,slab = data.T
    sx = x[0]
    sy = y[0]
    new_cord=np.zeros(len(x))
    for uu in range(1,len(x)):
        new_cord[uu]=fd.getDistance(y[uu], x[uu], sy, sx)
    ax3.plot(new_cord,slab/1000,color='green',lw=5,label='observation')
    ax3.legend(fontsize=24)
    fig2.savefig(figpath+'gravity_topo_phase'+str(model)+'_'+str(frame)+'.png')
