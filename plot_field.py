#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 20:59:46 2022

@author: chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
#import gravity as fg
import matplotlib
import matplotlib as mpl
#matplotlib.use('Agg')
from matplotlib import cm
from netCDF4 import Dataset
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
#import flac_interpolate as fi


plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- DO WHAT -----------------------------------
# Model
Cocos           = 0
Nazca           = 1
### pdf or png
png             = 1
pdf             = 0

### plot
shot_12         = 0
shot_5          = 0  ### phase, viscosity, density, dynamic pressure, srII
plot_phase      = 0
plot_density    = 0
plot_pressure   = 1
plot_sxx        = 0
plot_phase_topo = 0
shot_3          = 0 ### phase, viscosity, dynamic pressure


#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Nazca/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'

if Cocos:
    xmin,xmax=450,1000
    zmin,zmax=-150,10
    model = 'Ref_Cocos'
# elif Nazca:
#     xmin,xmax=250,1000
#     zmin,zmax=-200,10
#     model = 'Ref_Nazca'
elif Nazca:
    xmin,xmax=250,1000
    zmin,zmax=-200,10
    model = 'Nazca_a0702'
frame = 55
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

    
g=10
#------------------------------------------------------------------------------
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
def get_vis(frame):
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    vis = fl.read_visc(frame)
    xtop,ztop = fd.get_topo(x,z)
    return x,z,ele_x,ele_z,vis,ztop
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre
#------------------------------------------------------------------------------
end=250
# listmin=np.zeros(end)
# listmax=np.zeros(end)
# for frame in range(1,150+1):
if shot_12:
    fig, (ax)= plt.subplots(4,3,figsize=(39,18))
    filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
    x,z,ph=np.loadtxt(filepath).T
    ax[0,0].pcolormesh(x,z,ph,cmap = phase15,vmax=20,vmin=1,shading = 'auto')
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    cx=ax[0,0].contour(x,z,temp,cmap = 'rainbow',levels =[200,400,600,800,1000,1200],linewidths=1)
    ax[0,0].contour(x,z,temp,levels =[1300],linewidths=2,colors = '#F08080',linestyles='dashed')
    ax[0,0].clabel(cx, inline=True, fontsize=10,colors='white',fmt="%1.0f")
    
    
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    
    cc = plt.cm.get_cmap('jet')
    cbvis=ax[1,0].pcolormesh(x,z,vis,cmap=cc,vmin=20, vmax=27)
    ax[1,0].set_title('viscosity',fontsize=25)
    cax = plt.axes([0.36, 0.518, 0.008, 0.165])
    cc1=fig.colorbar(cbvis, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
   
    
    den = fl.read_density(frame)
    cc = plt.cm.get_cmap('RdBu_r')
    cbden=ax[2,0].pcolormesh(x,z,den,cmap=cc,vmin=2800,vmax=3400)
    ax[2,0].set_title('density',fontsize=25)
    cax = plt.axes([0.36, 0.32, 0.008, 0.165])
    cc1=fig.colorbar(cbden, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    srii = fl.read_srII(frame)
    cbsrii = plt.cm.get_cmap('afmhot')
    cbsrii=ax[3,0].pcolormesh(x,z,srii,cmap=cbsrii,vmin=-15,vmax=-12)
    ax[3,0].set_title('srII',fontsize=25)
    cax = plt.axes([0.36, 0.126, 0.008, 0.165])
    cc1=fig.colorbar(cbsrii, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    ck=ax[0,1].pcolormesh(x,z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax[0,1].set_title('pressure',fontsize=25)
    cax = plt.axes([0.635, 0.708, 0.008, 0.165])
    cc1=fig.colorbar(ck, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    fmelt=fl.read_fmelt(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    ck=ax[1,1].pcolormesh(x,z,fmelt*100,cmap=ck,vmin=0,vmax=3)
    ax[1,1].set_title('fmelt',fontsize=25)
    cax = plt.axes([0.635, 0.518, 0.008, 0.165])
    cc1=fig.colorbar(ck, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    fmagma = fl.read_fmagma(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    ck=ax[2,1].pcolormesh(x,z,fmagma*100,cmap=ck,vmin=0,vmax=0.03)
    ax[2,1].set_title('fmagma',fontsize=25)
    cax = plt.axes([0.635, 0.32, 0.008, 0.165])
    cc1=fig.colorbar(ck, ax=ax[3],cax=cax)
    cc1.ax.tick_params(labelsize=20)


    sxx = fl.read_sxx(frame)
    cbsxx = plt.cm.get_cmap('afmhot')
    cbsxx=ax[3,1].pcolormesh(x,z,sxx,cmap=cbsxx,vmin=-5,vmax=5)
    ax[3,1].set_title('sxx',fontsize=25)
    cax = plt.axes([0.635, 0.126, 0.008, 0.165])
    cc1=fig.colorbar(cbsxx, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    sii = fl.read_sII(frame)
    cbsii = plt.cm.get_cmap('afmhot_r')
    cbsii=ax[0,2].pcolormesh(x,z,sii,cmap=cbsii,vmin=0,vmax=9)
    ax[0,2].set_title('sII',fontsize=25)
    cax = plt.axes([0.91, 0.708, 0.008, 0.165])
    cc1=fig.colorbar(cbsii, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    sxz = fl.read_sxz(frame)
    cbsxz = plt.cm.get_cmap('afmhot')
    cbsxz=ax[1,2].pcolormesh(x,z,sxz,cmap=cbsxz,vmin=-4,vmax=4)
    ax[1,2].set_title('sxz',fontsize=25)
    cax = plt.axes([0.91, 0.518, 0.008, 0.165])
    cc1=fig.colorbar(cbsxz, ax=ax[1],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    szz=fl.read_szz(frame)
    cbszz = plt.cm.get_cmap('afmhot')
    cbszz=ax[2,2].pcolormesh(x,z,szz,cmap=cbszz,vmin=-4,vmax=4)
    ax[2,2].set_title('szz',fontsize=25)
    cax = plt.axes([0.91, 0.32, 0.008, 0.165])
    cc1=fig.colorbar(cbszz, ax=ax[3],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    
    
    strain=fl.read_strain(frame)
    cbszz = plt.cm.get_cmap('afmhot')
    cbszz=ax[3,2].pcolormesh(x,z,strain[0],cmap=cbszz,vmin=-2,vmax=2)
    ax[3,2].set_title('strain',fontsize=25)
    cax = plt.axes([0.91, 0.126, 0.008, 0.165])
    cc1=fig.colorbar(cbszz, ax=ax[3],cax=cax)
    cc1.ax.tick_params(labelsize=20)
   
    
    for yy in range(len(ax)):
        for qq in range(len(ax[0])):
            ax[yy,qq].set_xlim(xmin,xmax)
            ax[yy,qq].set_ylim(zmin,zmax)
            # ax[0,0].set_title(str(model)+' at '+str(round(fl.time[frame-1],1))+' Myr',fontsize=24)
            ax[yy,0].set_ylabel('Depth (km)',fontsize=20)
            ax[3,qq].set_xlabel('Distance (km)',fontsize=20)
            ax[yy,qq].set_aspect('equal')
            bwith = 3
            ax[yy,qq].spines['bottom'].set_linewidth(bwith)
            ax[yy,qq].spines['top'].set_linewidth(bwith)
            ax[yy,qq].spines['right'].set_linewidth(bwith)
            ax[yy,qq].spines['left'].set_linewidth(bwith)
            # ax[yy,0].spines['top'].set_visible(False)
            ax[yy,qq].tick_params(axis='x', labelsize=16 )
            ax[yy,qq].tick_params(axis='y', labelsize=16 )
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_allfield_150.png')
    # if pdf:
    #     if Cocos:
    #         fig.savefig(figpath+'Ref_Cocos/'+model+'frame_'+str(frame)+'_snapshot.pdf')
    #     if Nazca:
    #         fig.savefig(figpath+'Ref_Nazca/'+model+'frame_'+str(frame)+'_snapshot.pdf')

if shot_5:
# import time
# frame_list=[26,51,76,101,126,150]
# frame_list=[150]
# for frame in range(1,end+1,5):
# for frame in frame_list:
    fig, (ax)= plt.subplots(5,1,figsize=(15,22))

    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax[0].pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    cx=ax[0].contour(x,-z,temp,cmap = 'rainbow',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax[0].contour(x,-z,temp,levels =[1300],linewidths=4,colors = '#F08080',linestyles='dashed')
    
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax[1].pcolormesh(x,-z,vis,cmap=cc,vmin=20, vmax=25)
    ax[1].contour(x,-z,temp,levels =[400,600,800],linewidths=4,colors = '#C0C0C0')
    cax = plt.axes([0.945, 0.60, 0.01, 0.121])
    cc1=fig.colorbar(cbvis, ax=ax[1],cax=cax)
    cc1.set_label(label='Log$_{10}$ Viscosity (Pa $\cdot$ s)', size=23)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    
    
    den = fl.read_density(frame)
    cc = plt.cm.get_cmap('RdBu_r')
    cbden=ax[2].pcolormesh(x,-z,den,cmap=cc,vmin=2800,vmax=3400)
    ax[2].contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cax = plt.axes([0.945, 0.442, 0.01, 0.121])
    cc2=fig.colorbar(cbden, ax=ax[2],cax=cax)
    cc2.set_label(label='Density (kg/m$^3$)', size=25)
    cc2.ax.tick_params(labelsize=20)
    cc2.ax.yaxis.set_label_position('left')
    melt=fl.read_fmelt(frame)
    cbb = plt.cm.get_cmap('spring_r')
    # cbmag=ax[2].contour(ele_x,-ele_z,melt*100,cmap=cbb,vmin=0.1,vmax=10)
    ax[2].scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='yellow')
    
    
        
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    cbpre=ax[3].pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax[3].contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cax = plt.axes([0.945, 0.285, 0.01, 0.121])
    cc3=fig.colorbar(cbpre, ax=ax[1],cax=cax)
    cc3.set_label(label='Pressure (MPa)', size=25)
    cc3.ax.tick_params(labelsize=20)
    cc3.ax.yaxis.set_label_position('left')
    
    srii = fl.read_srII(frame)
    cbsrii = plt.cm.get_cmap('magma_r')
    cbsrii=ax[4].pcolormesh(x,-z,srii,cmap=cbsrii,vmin=-15,vmax=-12)
    ax[4].contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cax = plt.axes([0.945, 0.13, 0.01, 0.121])
    cc1=fig.colorbar(cbsrii, ax=ax[1],cax=cax)
    cc1.set_label(label='srII (s$^{-1}$)', size=25)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    

    ax[-1].set_xlabel('Distance (km)',fontsize=25)
    ax[0].set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=20)
    for yy in range(len(ax)):
        ax[yy].set_xlim(xmin,xmax)
        ax[yy].set_ylim(-zmin,-zmax)
        ax[yy].set_ylabel('Depth (km)',fontsize=25)
        
        ax[yy].set_aspect('equal')
        bwith = 3
        ax[yy].spines['bottom'].set_linewidth(bwith)
        ax[yy].spines['top'].set_linewidth(bwith)
        ax[yy].spines['right'].set_linewidth(bwith)
        ax[yy].spines['left'].set_linewidth(bwith)
        # ax[yy,0].spines['top'].set_visible(False)
        ax[yy].tick_params(axis='x', labelsize=22)
        ax[yy].tick_params(axis='y', labelsize=22)
    
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_5field_'+str(-zmin)+'.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_5field_'+str(-zmin)+'.pdf')
    # fig.gca()
    # plt.close(fig)
   
    
    
# for frame in range(2,end+1):
if plot_phase:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    from netCDF4 import Dataset
    # data = Dataset(savepath+model+'_phase3.'+str(frame)+'.grd', mode='r')
    # data = Dataset(savepath+model+'.grd', mode='r')
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax.pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(200,0,num=5)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylabel('Depth (km)',fontsize=30)
    ax.set_xlabel('Distance (km)',fontsize=30)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    #ax.set_xlim(250,1000)
    xmajor_ticks = np.linspace(250,1000,num=7)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=30)
    fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    # fig.gca()
    # plt.close(fig)
if plot_density:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    den = fl.read_density(frame)
    cc = plt.cm.get_cmap('seismic_r')
    cbden=ax.pcolormesh(x,-z,den,cmap=cc,vmin=2900,vmax=3500)
    cax = plt.axes([0.945, 0.382, 0.01, 0.251])
    cc2=fig.colorbar(cbden, ax=ax,cax=cax)
    cc2.set_label(label='Density (kg/m$^3$)', size=25)
    cc2.ax.tick_params(labelsize=20)
    cc2.ax.yaxis.set_label_position('left')
    melt=fl.read_fmelt(frame)
    cbb = plt.cm.get_cmap('spring_r')
    # cbmag=ax.contour(ele_x,-ele_z,melt*100,cmap=cbb,vmin=0.1,vmax=10)
    ax.contour(xt,-zt,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    xmajor_ticks = np.linspace(250,1000,num=7)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
frame =  10*5+1   
#if plot_pressure:
for frame in range(1,250+1):
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    cbpre=ax.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.945, 0.385, 0.01, 0.231])
    cc3=fig.colorbar(cbpre, ax=ax,cax=cax)
    cc3.set_label(label='Pressure (MPa)', size=25)
    cc3.ax.tick_params(labelsize=20)
    cc3.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    #ax.set_xlim(250,1000)
    ax.set_ylim(-zmin,-zmax)
    xmajor_ticks = np.linspace(250,1000,num=7)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    # ax.set_xlim(xmin,750)
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=30)
    fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
if plot_sxx:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    sxx = fl.read_sxx(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
    # cbsxx=ax.pcolormesh(x,-z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='flat')
    # cbsxx=ax.pcolormesh(x,-z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='auto')
    ax.set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_xlim(xmin,xmax)
    xmajor_ticks = np.linspace(250,1000,num=7)
    ax.set_xticks(xmajor_ticks)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    fig.savefig(figpath+model+'frame_'+str(frame)+'_sxx_gourand.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
if plot_phase_topo:
# for frame in range(2,229):
    fig, (ax0,ax1)= plt.subplots(2,1,figsize=(17,8),gridspec_kw={'height_ratios':[2,8]})
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- Topography -----------------------------
    xtop,ztop=fd.get_topo(xt,zt)
    ref_top = ztop[-12]
    ax0.plot(xtop,ztop,lw=2)
    ax0.axhline(y=-1,xmin=0,xmax=1,color='0.6',lw=2,ls='--')
    #--------------------- phase plotting -------------------------
    from netCDF4 import Dataset
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax1.pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    ax1.contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    for ax in (ax0,ax1):
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.tick_params(axis='x', labelsize=23)
        ax.tick_params(axis='y', labelsize=23)
        ax.set_xlim(250,1000)
        ax.set_xlim(xmin,xmax)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax.set_xticks(xmajor_ticks)
    ax0.set_ylim(-10,3)
    ymajor_ticks = np.linspace(200,0,num=5)
    ax1.set_yticks(ymajor_ticks)
    ax1.set_ylabel('Depth (km)',fontsize=30)
    ax1.set_xlabel('Distance (km)',fontsize=30)
    ax1.set_ylim(200,-30)
    ax1.set_ylim(-zmin,-zmax)
    
    
    ax1.set_aspect('equal')
    ax0.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=30)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    # fig.gca()
    # plt.close(fig)
if shot_3:
#for frame in range(200,250+1):
    fig, (ax)= plt.subplots(3,1,figsize=(15,12))
    
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax[0].pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
    cx=ax[0].contour(x,-z,temp,cmap = 'rainbow',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax[0].contour(x,-z,temp,levels =[1300],linewidths=4,colors = '#F08080',linestyles='dashed')
    melt=fl.read_fmelt(frame)
    # cbb = plt.cm.get_cmap('spring_r')
    # cbmag=ax[0].contour(ele_x,-ele_z,melt*100,colors='w',vmin=0.1,vmax=10)
    ax[0].scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e3*3,c='w')
    
    
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax[1].pcolormesh(x,-z,vis,cmap=cc,vmin=20, vmax=24)
    ax[1].contour(x,-z,temp,levels =[400,600,800],linewidths=4,colors = '#C0C0C0')
    cax = plt.axes([0.875, 0.38, 0.015, 0.231])
    cc1=fig.colorbar(cbvis, ax=ax[1],cax=cax)
    cc1.set_label(label='Log$_{10}$ Viscosity (Pa $\cdot$ s)', size=20)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    
        
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    cbpre=ax[2].pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax[2].contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cax = plt.axes([0.875, 0.125, 0.015, 0.221])
    cc3=fig.colorbar(cbpre, ax=ax[1],cax=cax)
    cc3.set_label(label='Pressure (MPa)', size=23)
    cc3.ax.tick_params(labelsize=20)
    cc3.ax.yaxis.set_label_position('left')
    
    
    ax[-1].set_xlabel('Distance (km)',fontsize=25)
    ax[0].set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    for yy in range(len(ax)):
        ax[yy].set_xlim(xmin,xmax)
        ax[yy].set_ylim(-zmin,-zmax)
        ax[yy].set_ylabel('Depth (km)',fontsize=25)
        
        ax[yy].set_aspect('equal')
        bwith = 3
        ax[yy].spines['bottom'].set_linewidth(bwith)
        ax[yy].spines['top'].set_linewidth(bwith)
        ax[yy].spines['right'].set_linewidth(bwith)
        ax[yy].spines['left'].set_linewidth(bwith)
        # ax[yy,0].spines['top'].set_visible(False)
        ax[yy].tick_params(axis='x', labelsize=22)
        ax[yy].tick_params(axis='y', labelsize=22)
    
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_3field_'+str(-zmin)+'.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_3field_'+str(-zmin)+'.pdf')
        fig.close()
