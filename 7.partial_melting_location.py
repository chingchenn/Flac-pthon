#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:19:28 2021

@author: ji-chingchen
"""
import os,sys
import flac
import numpy as np
import pandas as pd
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
from numpy import unravel_index
import matplotlib.pyplot as plt

fig1 = 0 ## single_bar_plot_melting
fig2 = 0 ## meltloc distance from trench
fig3 = 0
fig4 = 1 

model=sys.argv[1]
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/04model/'
#path = 'F:/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
bwith = 3

def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
def melting_location(start_vts=1,model_steps=end-1):
    melt=np.zeros(end)
    x_melt=np.zeros(end)
    z_melt= np.zeros(end)
    for i in range(1,end):
        x,z=fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        melt[i] = np.max(mm)
        maxindex_x=unravel_index(mm.argmax(),mm.shape)[0]
        maxindex_z=unravel_index(mm.argmax(),mm.shape)[1]
        ele_x, ele_z = nodes_to_elements(x,z)
        tren_x = (ele_z[:,0]).argmin()
        x_melt[i] = ele_x[maxindex_x,0]-ele_x[tren_x,0]
        z_melt[i]=fd.read_depth(z,maxindex_x,maxindex_z)
    return melt,x_melt,z_melt,time
def melting_phase():
    melt_num = np.zeros(end)
    phase_p3=np.zeros(end)
    phase_p4=np.zeros(end)
    phase_p9=np.zeros(end)
    phase_p10=np.zeros(end)
    for i in range(1,end):
        c=0;p9=0;p4=0;p10=0;p3=0
        mm=fl.read_fmelt(i)
        phase=fl.read_phase(i)
        area = fl.read_area(i)
        for xx in range(len(mm)):
            for zz in range(len(mm[0])):
                if mm[xx,zz] != 0:
                    if phase[xx,zz]==9:
                        p9 += area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==4:
                        p4 +=area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==10 or phase[xx,zz]==5:
                        p10 += area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==3:
                        p3 += area[xx,zz]*mm[xx,zz]/1e6
                    c +=1
        pk=c-p4-p9-p10-p3
        melt_num[i]=c
        phase_p3[i]=p3
        phase_p4[i]=p4
        phase_p9[i]=p9
        phase_p10[i]=p10
    return fl.time,phase_p3,phase_p4,phase_p9,phase_p10
if fig1:
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p9,phase_p10=melting_phase()
    fs.save_5txt(name,savepath,time,phase_p3,phase_p4,phase_p9,phase_p10)
    fig, (ax) = plt.subplots(1,1,figsize=(12,4))
    time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(path+name+'.txt').T
    ax.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
    ax.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments+basalt')
    ax.set_xlim(0,30)
    #ax.grid()
    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.legend(fontsize=25)
    #ax.set_title('Model : '+model,fontsize=25)
    ax.set_xlabel('Time (Myr)',fontsize=26)
    ax.set_ylabel('molten rocks (km3/km)',fontsize=26)
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    figpath = '/home/jiching/geoflac/figure/'
    fig.savefig(figpath+model+'_single_bar_plot_melting.png')
if fig2:
    name='metloc_for_'+model
    melt,xmelt,zmelt,time=melting_location()
    fs.save_3txt(name,savepath,time,melt,xmelt)
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    qqq=ax.scatter(time[melt>0],xmelt[melt>0],c=melt[melt>0],cmap='OrRd',vmin=0.0,vmax=0.05)
    cbar=fig.colorbar(qqq,ax=ax)
    ax.set_ylim(0,400)
    ax.set_xlim(0,30)
    ax.set_title(str(model)+" Melting location",fontsize=24)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    cbar.set_label('Melting %',fontsize=20)
    ax.set_ylabel('Distance with trench (km)',fontsize=20)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_metloc.png')
if fig3:
    melt,xmelt,zmelt,time=melting_location()
    fs.save_3txt('metdep_for_'+model,savepath,time,melt,zmelt)
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    time,melt,zmelt=np.loadtxt(savepath+'metdep_for_'+model+'.txt').T
    qqq=ax.scatter(time[melt>0],zmelt[melt>0],c=melt[melt>0],cmap='OrRd',vmin=0.0,vmax=0.05)
    cbar=fig.colorbar(qqq,ax=ax)
    ax.set_ylim(150,0)
    ax.set_xlim(0,30)
    ax.set_title(str(model)+" Melting location",fontsize=24)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    cbar.set_label('Melting %',fontsize=20)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_metdep.png')
if fig4:
    rainbow = cm.get_cmap('gray_r',end)
    meltcolor = cm.get_cmap('OrRd',end)
    newcolors = rainbow(np.linspace(0, 1, end))
    time_color = meltcolor(np.linspace(0,1,end))
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = nodes_to_elements(x,z)
        magma_chamber = fl.read_fmagma(i)
        melt = fl.read_fmelt(i)
        ax.scatter(ele_x[magma_chamber>1e-4],-ele_z[magma_chamber>1e-4],color=newcolors[i],zorder=1,s=10)
        if len(ele_x[melt>1e-4]) !=0:
            time = fl.time[i]
            qqq=ax.scatter(ele_x[melt>1e-4],-ele_z[melt>1e-4],color=time_color[i],s = 10)
    ax.set_ylim(150,0)
    #ax.set_xlim(0,1200)
    ax.set_title(str(model)+" Melting location",fontsize=24)
    ax.set_xlabel('X location',fontsize=20)
    ax.set_ylabel('Depth (km)',fontsize=20)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_melting_location_2D.png')
