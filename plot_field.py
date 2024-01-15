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
### make phase.grd file in 34.phase_interpolate.py

plt.rcParams["font.family"] = "Arial"
#---------------------------------- DO WHAT -----------------------------------
# Model
Cocos           = 0
Nazca           = 1
### pdf or png
png             = 0
pdf             = 0

### plot
shot_12         = 0
shot_5          = 0  ### phase, viscosity, density, dynamic pressure, srII
plot_phase      = 0
plot_viscosity  = 0
plot_density    = 0
plot_pressure   = 0
plot_sxx        = 1
plot_sII        = 0
plot_srII       = 0
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
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
# figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/AGU/POSTER/Poster_figure/'

if Cocos:
    xmin,xmax=500,900
    zmin,zmax=-120,10
    model = 'Ref_Cocos'
    shift = 550
    # model = 'Cocos11'
    # xmin,xmax=400,900
    # zmin,zmax=-150,10
# elif Nazca:
#     xmin,xmax=250,1000
#     zmin,zmax=-200,10
#     model = 'Ref_Nazca'
elif Nazca:
    xmin,xmax=250,1000
    zmin,zmax=-200,10
    model = 'Nazca_aa06'
    shift = 320
    #model = 'Ref_Nazca'
frame = 51
frame = 11
frame = 190
#frame = 200
# frame = 151
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

    
g=10
#------------------------------------------------------------------------------
def plot_snapshot(frame):
    x,z = fl.read_mesh(frame)
    xtop,ztop = fd.get_topo(x,z)
    phase = fl.read_phase(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    temp = fl.read_temperature(frame)
    return x, z, ele_x, ele_z, phase, temp, ztop
def get_vis(frame):
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    vis = fl.read_visc(frame)
    xtop,ztop = fd.get_topo(x,z)
    return x,z,ele_x,ele_z,vis,ztop
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
colors2=[
 '#C98F49', '#92C0DF', '#2553C7', '#FFFFFF', '#6495ED',
 '#2E8B57', '#524B52', '#9A32CD', '#6B8E23','#D4DBF4',
 '#D8BFD8','#999999','#F2C85B','#92C0DF','#999999',
 '#4CC552','#999999','#999999','#999999','#999999']
phase8= matplotlib.colors.ListedColormap(colors2)

def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre

def compute_s1(sxx, szz, sxz):
    mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
    theta = 0.5 * np.arctan2(2*sxz,  sxx-szz)

    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    nx, nz = sxx.shape
    tmp = np.zeros((nx, nz, 3), dtype=sxx.dtype)
    tmp[:,:,0] = mag * np.sin(theta)
    tmp[:,:,1] = mag * np.cos(theta)
    return tmp
#------------------------------------------------------------------------------
end=250
# listmin=np.zeros(end)
# listmax=np.zeros(end)
# for frame in range(1,150+1):
if shot_5:
# import time
# frame_list=[26,51,76,101,126,150]
# frame_list=[150]
# for frame in range(1,end+1,5):
# for frame in frame_list:
    fig, (ax)= plt.subplots(5,1,figsize=(15,22))

    file=model+'_phase_'+str(frame)+'.grd'
    # file=model+'_frame'+str(frame)+'_phase.grd'
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
   
   
    

if plot_phase:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    from netCDF4 import Dataset
    #file=model+'_frame'+str(frame)+'_phase.grd'
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=2)
    #--------------------- melting plotting -------------------------
    melt=fl.read_fmelt(frame)
    ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='w')
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(200,0,num=5)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylabel('Depth (km)',fontsize=26)
    # ax.set_xlabel('Distance (km)',fontsize=30)
    if Nazca:
        xmajor_ticks = np.linspace(250,1000,num=7)
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    ax.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    # fig.gca()
    # plt.close(fig)
    
if plot_viscosity:
# frame_list=[26,51,76,101,126,150]
# frame_list=[150]
# for frame in range(1,end+1,5):
# for frame in frame_list:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=25,shading='gouraud')
    cax = plt.axes([0.945, 0.382, 0.01, 0.251])
    cc1=fig.colorbar(cbvis, ax=ax,cax=cax)
    cc1.set_label(label='Log$_{10}$ Viscosity (Pa $\cdot$ s)', size=23)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    
    # cc2=fig.colorbar(cbden, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(xt,-zt,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.set_ylabel('Depth (km)',fontsize=26)
    # ax.set_xlabel('Distance (km)',fontsize=26)
    ax.tick_params(labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    if Nazca:
        xmajor_ticks = np.linspace(250,1000,num=7)
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_viscosity.png')
    if pdf:    
        fig.savefig(figpath+model+'frame_'+str(frame)+'_viscosity.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
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
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    if Nazca:
        xmajor_ticks = np.linspace(250,1000,num=7)
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-zmin,-zmax)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
#frame =  10*5+1   
if plot_pressure:
# for frame in [11,51,76,101,126,151]:
# for frame in range(1,250+1):
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
    ax.set_ylabel('Depth (km)',fontsize=26)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(labelsize=25)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    #xmajor_ticks = np.linspace(250,1000,num=6)
    #ax.set_xticks(xmajor_ticks)
    #ax.set_xlim(250,1000)
    ax.set_ylim(-zmin,-zmax)
    if Nazca:
        xmajor_ticks = np.linspace(250,1000,num=7)
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xlim(xmin,xmax)
    # ax.set_xlim(xmin,750)
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
if plot_sxx:
    skip = (slice(None, None, 5), slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx*100,cmap=cbsxx,vmin=-2000,vmax=2000,shading='gouraud')
    
    ax.quiver(ele_x[skip],-ele_z[skip],s1.T[skip],s3.T[skip])#,
             #angles='xy', scale_units='xy', scale=0.1,headwidth=3)
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
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(230,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(300,950)
    #xmajor_ticks = np.linspace(250,1000,num=7)
    #ax.set_xticks(xmajor_ticks)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=5)
       
    xxm,zzm,xtm = np.loadtxt(savepath+'Nazca_aa06_40_final_moho_slab.txt').T
    xxm=xxm[zzm<0]
    zzm=zzm[zzm<0]
    zzm = fd.moving_window_smooth(zzm,3)
    ax.plot(xxm+shift-5,-zzm+2,color='k',lw=5)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    #fig.savefig(figpath+model+'frame_'+str(frame)+'_sxx_gourand.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
if plot_sII:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sII = fl.read_sII(frame)*100
    cbsxx = plt.cm.get_cmap('Greens')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sII,cmap=cbsxx,vmin=0,vmax=600,shading='gouraud')
    ax.set_title('sII',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='sII (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='#A52A2A',levels =[200,400,600,800,1000,1200],linewidths=2)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_xlim(xmin,xmax)
    # xmajor_ticks = np.linspace(250,1000,num=7)
    # ax.set_xticks(xmajor_ticks)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=5)
    
    #fig.savefig(figpath+model+'frame_'+str(frame)+'_sII_gourand.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))

if plot_srII:
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    srII = fl.read_srII(frame)
    cbsxx = plt.cm.get_cmap('Blues')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,srII,cmap=cbsxx,vmin=-14,vmax=-12.5,shading='gouraud')
    #ax.scatter(ele_x[srII>-13.2],-ele_z[srII>-13.2],color='b')
    sII=fl.read_sII(frame)*100
    print(np.average(sII[srII>-13.2]))
    ax.set_title('srII',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='srII', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='#A52A2A',levels =[200,400,600,800,1000,1200],linewidths=2)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_xlim(xmin,xmax)
    # xmajor_ticks = np.linspace(250,1000,num=7)
    # ax.set_xticks(xmajor_ticks)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    #xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=5)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    #fig.savefig(figpath+model+'frame_'+str(frame)+'_sII_gourand.pdf')
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
        ax.tick_params(axis='x', labelsize=23)
        ax.tick_params(axis='y', labelsize=23)
        ax.set_xlim(xmin,xmax)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax.set_xticks(xmajor_ticks)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(bwith)
    ax0.set_ylim(-10,3)
    ymajor_ticks = np.linspace(200,0,num=5)
    ax1.set_yticks(ymajor_ticks)
    ax1.set_ylabel('Depth (km)',fontsize=26)
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
# for frame in range(2,250+1):
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
        ax[yy].set_ylabel('Depth (km)',fontsize=26)
        ax[yy].set_aspect('equal')
        bwith = 3
        # ax[yy,0].spines['top'].set_visible(False)
        ax[yy].tick_params(axis='x', labelsize=22)
        ax[yy].tick_params(axis='y', labelsize=22)
        for axis in ['top','bottom','left','right']:
            ax[yy].spines[axis].set_linewidth(bwith)
    
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_3field_'+str(-zmin)+'.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_snapshot_3field_'+str(-zmin)+'.pdf')
        fig.close()
