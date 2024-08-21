#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 16:54:09 2022

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
png             = 0
pdf             = 0

### plot
plot_sxx        = 0
plot_srII       = 1
plot_sxz        = 0
plot_pressure   = 0
plot_phase      = 0

plot_Nazca      = 0
plot_Cocos      = 0
plot_Nazca_h    = 0
plot_Cocos_h    = 0


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
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/AGU/POSTER/Poster_figure/'

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
    
g=10
fig, (ax)= plt.subplots(2,1,figsize=(16,14),gridspec_kw={'height_ratios':[3,2]})
frame = 151
bwith = 3

def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre

if plot_sxx:
    #--------------------- plotting -------------------------
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    sxx = fl.read_sxx(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[0].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-200,vmax=200,shading='gouraud')
    ax[0].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    #--------------------- plotting -------------------------
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    sxx = fl.read_sxx(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[1].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-200,vmax=200,shading='gouraud')
    ax[1].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    for qq in ax:
        qq.set_aspect('equal')
        qq.spines['bottom'].set_linewidth(bwith)
        qq.spines['top'].set_linewidth(bwith)
        qq.spines['right'].set_linewidth(bwith)
        qq.spines['left'].set_linewidth(bwith)
        qq.tick_params(axis='x', labelsize=23)
        qq.tick_params(axis='y', labelsize=23)
        qq.set_ylabel('Depth (km)',fontsize=28)
        qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
    xmajor_ticks = np.linspace(300,900,num=13)
    ax[0].set_xticks(xmajor_ticks)
    ax[0].set_ylim(150,-10)
    ax[0].set_xlim(300,900)
    ax[1].set_ylim(100,-10)
    ax[1].set_xlim(500,850)
    ax[-1].set_xlabel('Distance (km)',fontsize=30)
    ax[0].set_title('Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    ax[1].set_title('Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    fig.savefig(figpath+'Sxx_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
if plot_srII:
    #--------------------- plotting -------------------------
    model = 'Nazca_aa06'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    srII = fl.read_srII(frame)
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[0].pcolormesh(ele_x,-ele_z,srII,cmap=cbsxx,vmin=-20,vmax=-11,shading='gouraud')
    ax[0].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    #--------------------- plotting -------------------------
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    sxx = fl.read_sxx(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[1].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-200,vmax=200,shading='gouraud')
    ax[1].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    for qq in ax:
        qq.set_aspect('equal')
        qq.spines['bottom'].set_linewidth(bwith)
        qq.spines['top'].set_linewidth(bwith)
        qq.spines['right'].set_linewidth(bwith)
        qq.spines['left'].set_linewidth(bwith)
        qq.tick_params(axis='x', labelsize=23)
        qq.tick_params(axis='y', labelsize=23)
        qq.set_ylabel('Depth (km)',fontsize=28)
        qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
    xmajor_ticks = np.linspace(300,900,num=13)
    ax[0].set_xticks(xmajor_ticks)
    ax[0].set_ylim(150,-10)
    ax[0].set_xlim(300,900)
    ax[1].set_ylim(100,-10)
    ax[1].set_xlim(500,850)
    ax[-1].set_xlabel('Distance (km)',fontsize=30)
    ax[0].set_title('Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    ax[1].set_title('Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    fig.savefig(figpath+'Sxx_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))

if plot_sxz:
    #--------------------- plotting -------------------------
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    sxx = fl.read_sxz(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[0].pcolormesh(x,-z,sxx,cmap=cbsxx,vmin=-200,vmax=200)
    ax[0].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='$\sigma_{xz}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    #--------------------- plotting -------------------------
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    sxx = fl.read_sxz(frame)*100
    cbsxx = plt.cm.get_cmap('seismic')
    cbsxx=ax[1].pcolormesh(x,-z,sxx,cmap=cbsxx,vmin=-200,vmax=200)
    ax[1].set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    for qq in ax:
        qq.set_aspect('equal')
        qq.spines['bottom'].set_linewidth(bwith)
        qq.spines['top'].set_linewidth(bwith)
        qq.spines['right'].set_linewidth(bwith)
        qq.spines['left'].set_linewidth(bwith)
        qq.tick_params(axis='x', labelsize=23)
        qq.tick_params(axis='y', labelsize=23)
        qq.set_ylabel('Depth (km)',fontsize=28)
        qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
    xmajor_ticks = np.linspace(300,900,num=13)
    ax[0].set_xticks(xmajor_ticks)
    ax[0].set_ylim(200,-10)
    ax[0].set_xlim(300,900)
    ax[1].set_ylim(100,-10)
    ax[1].set_xlim(500,850)
    ax[-1].set_xlabel('Distance (km)',fontsize=30)
    ax[0].set_title('Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    ax[1].set_title('Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    # fig.savefig(figpath+'Sxz_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
    
if plot_pressure:
    #--------------------- plotting -------------------------
    model = 'Ref_Nazca'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    temp = fl.read_temperature(frame)
    x,z,new_pre = dynamics_pressure(frame)
    ck = plt.cm.get_cmap('RdYlBu_r')
    cbpre=ax[0].pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    cax = plt.axes([0.945, 0.315, 0.01, 0.371])
    cc1=fig.colorbar(cbpre, ax=ax[0],cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='Pressure (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    #--------------------- plotting -------------------------
    model = 'Cocos_a0646'
    os.chdir(path+model)
    fl = flac.Flac()
    temp = fl.read_temperature(frame)
    x,z = fl.read_mesh(frame)
    sxx = fl.read_sxz(frame)
    x,z,new_pre = dynamics_pressure(frame)
    cbsxx=ax[1].pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
    ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    for qq in ax:
        qq.set_aspect('equal')
        qq.spines['bottom'].set_linewidth(bwith)
        qq.spines['top'].set_linewidth(bwith)
        qq.spines['right'].set_linewidth(bwith)
        qq.spines['left'].set_linewidth(bwith)
        qq.tick_params(axis='x', labelsize=23)
        qq.tick_params(axis='y', labelsize=23)
        qq.set_ylabel('Depth (km)',fontsize=28)
        qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
    xmajor_ticks = np.linspace(300,900,num=13)
    ax[0].set_xticks(xmajor_ticks)
    ax[0].set_ylim(200,-10)
    ax[0].set_xlim(300,900)
    ax[1].set_ylim(100,-10)
    ax[1].set_xlim(500,850)
    ax[-1].set_xlabel('Distance (km)',fontsize=30)
    ax[0].set_title('Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    ax[1].set_title('Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    fig.savefig(figpath+'Pressure_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
if plot_phase:
    #--------------------- phase plotting -------------------------
    model = 'Ref_Nazca'
    os.chdir(path+model)
    fl = flac.Flac()
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    from netCDF4 import Dataset
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax[0].pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    ax[0].contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
    #--------------------- phase plotting -------------------------
    model = 'Cocos_a0646'
    os.chdir(path+model)
    fl = flac.Flac()
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    phh=ph.data[ph.data>0]
    ax[1].pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
    ax[1].contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
    
    # ---------------------- plot setting --------------------------
    for qq in ax:
        qq.set_aspect('equal')
        qq.spines['bottom'].set_linewidth(bwith)
        qq.spines['top'].set_linewidth(bwith)
        qq.spines['right'].set_linewidth(bwith)
        qq.spines['left'].set_linewidth(bwith)
        qq.tick_params(axis='x', labelsize=23)
        qq.tick_params(axis='y', labelsize=23)
        qq.set_ylabel('Depth (km)',fontsize=28)
        qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
    xmajor_ticks = np.linspace(300,900,num=13)
    ax[0].set_xticks(xmajor_ticks)
    ax[0].set_ylim(200,-10)
    ax[0].set_xlim(300,900)
    ax[1].set_ylim(100,-10)
    ax[1].set_xlim(500,850)
    
    ax[-1].set_xlabel('Distance (km)',fontsize=30)
    ax[0].set_title('Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    ax[1].set_title('Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
    # fig.savefig(figpath+'Phase_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
if plot_Nazca:
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    for frame in [101,126,151,176,201]:
        fig, (ax)= plt.subplots(2,1,figsize=(16,14),gridspec_kw={'height_ratios':[1,1]})
        temp = fl.read_temperature(frame)
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax[0].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax[0].set_title('sxx',fontsize=25)
        cax = plt.axes([0.945, 0.565, 0.01, 0.271])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        x,z,new_pre = dynamics_pressure(frame)
        ck = plt.cm.get_cmap('RdYlBu_r')
        cbpre=ax[1].pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
        cax = plt.axes([0.945, 0.155, 0.01, 0.271])
        cc1=fig.colorbar(cbpre, ax=ax[1],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='Pressure (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        # ---------------------- plot setting --------------------------
        xmajor_ticks = np.linspace(300,1000,num=8)
        for qq in ax:
            qq.set_aspect('equal')
            qq.spines['bottom'].set_linewidth(bwith)
            qq.spines['top'].set_linewidth(bwith)
            qq.spines['right'].set_linewidth(bwith)
            qq.spines['left'].set_linewidth(bwith)
            qq.tick_params(axis='x', labelsize=23)
            qq.tick_params(axis='y', labelsize=23)
            qq.set_ylabel('Depth (km)',fontsize=28)
            qq.set_ylim(200,-10)
            qq.set_xlim(300,900)
            qq.set_xticks(xmajor_ticks)
            qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
    
        ax[-1].set_xlabel('Distance (km)',fontsize=30)
        ax[0].set_title('$\sigma_{xx}$ of Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        ax[1].set_title('Dynamics pressure of Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        # fig.savefig(figpath+'Sxx_and_Pressure_Chile_of_'+str(frame)+'.pdf')
    
if plot_Cocos:
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    for frame in [101,126,151,176,201]:
        fig, (ax)= plt.subplots(3,1,figsize=(16,21),gridspec_kw={'height_ratios':[1,1,1]})
        temp = fl.read_temperature(frame)
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax[0].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax[0].set_title('sxx',fontsize=25)
        cax = plt.axes([0.945, 0.675, 0.01, 0.191])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        x,z,new_pre = dynamics_pressure(frame)
        ck = plt.cm.get_cmap('RdYlBu_r')
        cbpre=ax[1].pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
        cax = plt.axes([0.945, 0.409, 0.01, 0.191])
        cc1=fig.colorbar(cbpre, ax=ax[1],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='Pressure (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        
        
        # model = 'Cocos_a0646'
        # os.chdir(path+model)
        # fl = flac.Flac()
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)
        file=model+'_phase_'+str(frame)+'.grd'
        data = Dataset(savepath+file, mode='r')
        x = data.variables['x'][:]
        z = data.variables['y'][:]
        ph = data.variables['z'][:]
        phh=ph.data[ph.data>0]
        ax[2].pcolormesh(x,-z,ph,cmap=phase15,vmin=1, vmax=20)
        ax[2].contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
                # ---------------------- plot setting --------------------------
        # xmajor_ticks = np.linspace(500,0,num=13)
        for qq in ax:
            qq.set_aspect('equal')
            qq.spines['bottom'].set_linewidth(bwith)
            qq.spines['top'].set_linewidth(bwith)
            qq.spines['right'].set_linewidth(bwith)
            qq.spines['left'].set_linewidth(bwith)
            qq.tick_params(axis='x', labelsize=23)
            qq.tick_params(axis='y', labelsize=23)
            qq.set_ylabel('Depth (km)',fontsize=28)
            qq.set_ylim(100,-10)
            qq.set_xlim(500,850)
            # qq.set_xticks(xmajor_ticks)
            qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
    
        ax[-1].set_xlabel('Distance (km)',fontsize=30)
        ax[0].set_title('$\sigma_{xx}$ of Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        ax[1].set_title('Dynamics pressure at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        ax[2].set_title('Phase profile at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)

        # fig.savefig(figpath+'Sxx_and_Pressure_Phase_Mexico_of_'+str(frame)+'.pdf')
        # 
if plot_Nazca_h:
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    for frame in [76,101,126,151,176,201]:
        fig, (ax)= plt.subplots(1,2,figsize=(32,7))
        temp = fl.read_temperature(frame)
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax[0].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax[0].set_title('sxx',fontsize=25)
        cax = plt.axes([0.499, 0.23, 0.01, 0.55])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        x,z,new_pre = dynamics_pressure(frame)
        ck = plt.cm.get_cmap('RdYlBu_r')
        cbpre=ax[1].pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
        cax = plt.axes([0.925, 0.23, 0.01, 0.55])
        cc1=fig.colorbar(cbpre, ax=ax[1],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='Pressure (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        # ---------------------- plot setting --------------------------
        xmajor_ticks = np.linspace(300,1000,num=8)
        for qq in ax:
            qq.set_aspect('equal')
            qq.spines['bottom'].set_linewidth(bwith)
            qq.spines['top'].set_linewidth(bwith)
            qq.spines['right'].set_linewidth(bwith)
            qq.spines['left'].set_linewidth(bwith)
            qq.tick_params(axis='x', labelsize=23)
            qq.tick_params(axis='y', labelsize=23)
            
            qq.set_xlabel('Distance (km)',fontsize=30)
            qq.set_ylim(200,-10)
            qq.set_xlim(300,900)
            qq.set_xticks(xmajor_ticks)
            qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
    
        ax[0].set_ylabel('Depth (km)',fontsize=28)
        ax[0].set_title('$\sigma_{xx}$ of Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        ax[1].set_title('Dynamics pressure of Chile model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        # fig.savefig(figpath+'Sxx_and_Pressure_Chile_of_'+str(frame)+'.pdf')
    
if plot_Cocos_h:
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    for frame in [101,126,151,176,201]:
        fig, (ax)= plt.subplots(1,2,figsize=(32,7))
        temp = fl.read_temperature(frame)
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax[0].pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax[0].set_title('sxx',fontsize=25)
        cax = plt.axes([0.499, 0.23, 0.01, 0.55])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        x,z,new_pre = dynamics_pressure(frame)
        ck = plt.cm.get_cmap('RdYlBu_r')
        cbpre=ax[1].pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
        cax = plt.axes([0.925, 0.23, 0.01, 0.55])
        cc1=fig.colorbar(cbpre, ax=ax[1],cax=cax)
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='Pressure (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        
        # ---------------------- plot setting --------------------------
        # xmajor_ticks = np.linspace(500,0,num=13)
        for qq in ax:
            qq.set_aspect('equal')
            qq.spines['bottom'].set_linewidth(bwith)
            qq.spines['top'].set_linewidth(bwith)
            qq.spines['right'].set_linewidth(bwith)
            qq.spines['left'].set_linewidth(bwith)
            qq.tick_params(axis='x', labelsize=23)
            qq.tick_params(axis='y', labelsize=23)
            qq.set_xlabel('Distance (km)',fontsize=30)
            
            qq.set_ylim(100,-10)
            qq.set_xlim(500,850)
            # qq.set_xticks(xmajor_ticks)
            qq.set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
        ax[0].set_ylabel('Depth (km)',fontsize=28)
        ax[0].set_title('$\sigma_{xx}$ of Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        ax[1].set_title('Dynamics pressure of Mexico model at '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=30)
        # fig.savefig(figpath+'Sxx_and_Pressure_Mexico_of_'+str(frame)+'.pdf')