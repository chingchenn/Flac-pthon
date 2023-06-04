#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 08:45:40 2023

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
plot_phase      = 1
plot_viscosity  = 0
plot_density    = 0
plot_pressure   = 0
plot_pressure3  = 0
plot_phase_pressure = 0
plot_viscosity_pressure = 0
plot_sxx        = 0
plot_phase_topo = 0

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
# figpath = '/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'

if Cocos:
    xmin,xmax = 500,900
    zmin,zmax = -150,10
    model = 'Ref_Cocos'
    frame1 = 50
    frame2 = 70
    frame3 = 110
    frame4 = 190
elif Nazca:
    xmin,xmax = 250,1000
    zmin,zmax = -200,10
    model = 'Nazca_a0702'
    frame1 = 30
    frame2 = 60
    frame3 = 120
    frame4 = 140

os.chdir(path+model)
fl = flac.Flac()
time=fl.time
bwith = 4
fontsize=35
labelsize=30

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
#------------------------------------------------------------------------------
if plot_phase:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,26))
    if Nazca:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,22))

    #----------------------------- FIG1 -----------------------------
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    if Nazca:
        ax.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    #----------------------------- FIG2 -----------------------------
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax2.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax2.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    if Nazca:
        ax2.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax2.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    #----------------------------- FIG3 -----------------------------
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax3.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax3.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    if Nazca:
        ax3.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax3.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    #----------------------------- FIG4 -----------------------------
    frame = frame4
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax4.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax4.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    if Nazca:
        ax4.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax4.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    # ---------------------- plot setting --------------------------
    if Nazca:
        # xmajor_ticks = np.linspace(250,1000,num=16)
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3,ax4]:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
        aa.set_aspect('equal')
        aa.set_yticks(ymajor_ticks)
        aa.set_ylabel('Depth (km)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
        aa.set_xlim(xmin,xmax)
        aa.set_ylim(-zmin,-zmax)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax4.set_xlabel('Distance (km)',fontsize=fontsize)
    # ax.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
        
if plot_pressure3:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,26))
    if Nazca:
        # fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,22))
        fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))

    #----------------------------- FIG1 -----------------------------
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdBu_r')
    cbpre=ax.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.945, 0.385, 0.01, 0.231])
    cc3=fig.colorbar(cbpre, ax=ax,cax=cax)
    cc3.set_label(label='Pressure (MPa)', size=25)
    cc3.ax.tick_params(labelsize=20)
    cc3.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    #----------------------------- FIG2 -----------------------------
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    cbpre=ax2.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax2.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax2.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax2.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)    
    #----------------------------- FIG3 -----------------------------
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    cbpre=ax3.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax3.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax3.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax3.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)    
    # #----------------------------- FIG4 -----------------------------
    # frame = frame4
    # xt,zt = fl.read_mesh(frame)
    # ele_x,ele_z = flac.elem_coord(xt, zt)
    # temp = fl.read_temperature(frame)
    # #--------------------- phase plotting -------------------------
    # x,z,new_pre = dynamics_pressure(frame)
    # cbpre=ax4.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    # # cax = plt.axes([0.945, 0.385, 0.01, 0.231])
    # # cc3=fig.colorbar(cbpre, ax=ax4,cax=cax)
    # # cc3.set_label(label='Pressure (MPa)', size=25)
    # # cc3.ax4.tick_params(labelsize=20)
    # # cc3.ax4.yaxis.set_label_position('left')
    # ax4.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    # if Nazca:
    #     ax4.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    # if Cocos:
    #     ax4.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)            
        
    # ---------------------- plot setting --------------------------    
    if Nazca:
        # xmajor_ticks = np.linspace(250,1000,num=16)
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3]:#,ax4]:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
        aa.set_aspect('equal')
        aa.set_yticks(ymajor_ticks)
        aa.set_ylabel('Depth (km)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
        aa.set_xlim(xmin,xmax)
        aa.set_ylim(-zmin,-zmax)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax3.set_xlabel('Distance (km)',fontsize=fontsize)
    # ax.set_title(model,fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.pdf')

if plot_pressure:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,26))
    if Nazca:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,22))
        # fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))

    #----------------------------- FIG1 -----------------------------
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    ck = plt.cm.get_cmap('RdBu_r')
    cbpre=ax.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.945, 0.385, 0.01, 0.231])
    cc3=fig.colorbar(cbpre, ax=ax,cax=cax)
    cc3.set_label(label='Pressure (MPa)', size=25)
    cc3.ax.tick_params(labelsize=20)
    cc3.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    #----------------------------- FIG2 -----------------------------
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax2.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax2.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax2.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax2.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)    
    #----------------------------- FIG3 -----------------------------
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax3.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax3.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax3.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax3.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)    
    # #----------------------------- FIG4 -----------------------------
    frame = frame4
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    # #--------------------- phase plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax4.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax4.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax4.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax4.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)            
        
    # ---------------------- plot setting --------------------------    
    if Nazca:
        # xmajor_ticks = np.linspace(250,1000,num=16)
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3,ax4]:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
        aa.set_aspect('equal')
        aa.set_yticks(ymajor_ticks)
        aa.set_ylabel('Depth (km)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
        aa.set_xlim(xmin,xmax)
        aa.set_ylim(-zmin,-zmax)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax3.set_xlabel('Distance (km)',fontsize=fontsize)
    # ax.set_title(model,fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure.pdf')
   
if plot_phase_pressure:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,26))
    if Nazca:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,22))
        # fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))

    #----------------------------- FIG1 -----------------------------
    ax5 = ax[1] 
    ax = ax[0]
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    ck = plt.cm.get_cmap('RdBu_r')
    cbpre=ax5.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.928, 0.285, 0.01, 0.431])
    cc3=fig.colorbar(cbpre, ax=ax5,cax=cax)
    cc3.set_label(label='Dynamics Pressure (MPa)', size=35)
    cc3.ax.tick_params(labelsize=30)
    cc3.ax.yaxis.set_label_position('left')
    ax5.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)    
    
    #----------------------------- FIG2 -----------------------------
    ax6 = ax2[1] 
    ax2 = ax2[0]
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax2.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax2.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    # ax2.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],c='#DC143C')
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax6.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax6.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    
    if Nazca:
        ax2.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax2.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    #----------------------------- FIG3 -----------------------------
    ax7 = ax3[1] 
    ax3 = ax3[0]
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax3.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax3.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax7.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax7.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax3.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax3.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    #----------------------------- FIG4 -----------------------------
    ax8 = ax4[1] 
    ax4 = ax4[0]
    frame = frame4
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax4.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax4.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax8.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax8.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax4.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    if Cocos:
        ax4.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize)
    
    # ---------------------- plot setting --------------------------
    if Nazca:
        # xmajor_ticks = np.linspace(250,1000,num=16)
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3,ax4,ax5,ax6,ax7,ax8]:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
        aa.set_aspect('equal')
        aa.set_yticks(ymajor_ticks)
        aa.set_ylabel('Depth (km)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
        aa.set_xlim(xmin,xmax)
        aa.set_ylim(-zmin,-zmax)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax4.set_xlabel('Distance (km)',fontsize=fontsize)
    ax8.set_xlabel('Distance (km)',fontsize=fontsize)
    # ax.set_title(model,fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase_pressure.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase_pressure.pdf')
    
    
if plot_viscosity_pressure:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,26))
    if Nazca:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,22))
        # fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))
   
    ax5 = ax[1] 
    ax = ax[0]
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    #--------------------- viscosity plotting -------------------------
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
    cax = plt.axes([0.93, 0.582, 0.01, 0.251])
    cc1=fig.colorbar(cbvis, ax=ax,cax=cax)
    cc1.set_label(label='Log$_{10}$ Viscosity (Pa $\cdot$ s)', size=35)
    cc1.ax.tick_params(labelsize=30)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    ck = plt.cm.get_cmap('RdBu_r')
    cbpre=ax5.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.93, 0.185, 0.01, 0.251])
    cc3=fig.colorbar(cbpre, ax=ax5,cax=cax)
    cc3.set_label(label='Dynamics Pressure (MPa)', size=35)
    cc3.ax.tick_params(labelsize=30)
    cc3.ax.yaxis.set_label_position('left')
    ax5.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    if Cocos:
        ax.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')    
    
    #----------------------------- FIG2 -----------------------------
    ax6 = ax2[1] 
    ax2 = ax2[0]
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    ##--------------------- viscosity plotting -------------------------
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax2.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
    ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax6.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax6.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    
    if Nazca:
        ax2.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    if Cocos:
        ax2.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    
    #----------------------------- FIG3 -----------------------------
    ax7 = ax3[1] 
    ax3 = ax3[0]
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    ##--------------------- viscosity plotting -------------------------
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax3.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
    ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax7.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax7.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax3.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    if Cocos:
        ax3.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    
    #----------------------------- FIG4 -----------------------------
    ax8 = ax4[1] 
    ax4 = ax4[0]
    frame = frame4
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    ##--------------------- viscosity plotting -------------------------
    x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
    cc = plt.cm.get_cmap('jet')
    cbvis=ax4.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
    ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    #--------------------- pressure plotting -------------------------
    x,z,new_pre = dynamics_pressure(frame)
    new_pre = new_pre-np.median(new_pre)
    cbpre=ax8.pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax8.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
    if Nazca:
        ax4.text(270,170,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    if Cocos:
        ax4.text(520,120,str(np.round(fl.time[frame-1],0))+' Myr',fontsize=fontsize, color = 'white')
    
    # ---------------------- plot setting --------------------------
    if Nazca:
        # xmajor_ticks = np.linspace(250,1000,num=16)
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3,ax4,ax5,ax6,ax7,ax8]:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
        aa.set_aspect('equal')
        aa.set_yticks(ymajor_ticks)
        aa.set_ylabel('Depth (km)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
        aa.set_xlim(xmin,xmax)
        aa.set_ylim(-zmin,-zmax)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    ax4.set_xlabel('Distance (km)',fontsize=fontsize)
    ax8.set_xlabel('Distance (km)',fontsize=fontsize)
    # ax.set_title(model,fontsize=30)
    if png:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_viscosity_pressure.png')
    if pdf:
        fig.savefig(figpath+model+'frame_'+str(frame)+'_viscosity_pressure.pdf')