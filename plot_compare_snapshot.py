#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:32:41 2022

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
shot_12         = 0
shot_5          = 0  ### phase, viscosity, density, dynamic pressure, srII
plot_phase      = 0
plot_viscosity  = 1
plot_density    = 0
plot_pressure   = 0
plot_sxx        = 0
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
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Cocos/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
frame_list=[40,60,80,100,120]
frame_list=[140,160,180,200,220]
frame_list=[50,100,125,150]
if Cocos:
    xmin,xmax=450,1000
    zmin,zmax=-150,10
    model = 'Ref_Cocos'
# elif Nazca:
#     xmin,xmax=250,1000
#     zmin,zmax=-200,10
#     model = 'Ref_Nazca'
elif Nazca:
    xmin,xmax=250,1100
    zmin,zmax=-200,10
    model = 'Nazca_a0645'
    # model = 'Ref_Nazca'
frame = 55
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

def get_vis(frame):
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    vis = fl.read_visc(frame)
    xtop,ztop = fd.get_topo(x,z)
    return x,z,ele_x,ele_z,vis,ztop
def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,ele_x,ele_z,dypre
g=10
bwith = 3
if plot_viscosity:
    # fig, (ax)= plt.subplots(3,1,figsize=(13,13))
    fig, (ax)= plt.subplots(len(frame_list),1,figsize=(13,18))
    cc = plt.cm.get_cmap('jet')
    #--------------------- viscosity plotting -------------------------
    for qq,frame in enumerate(frame_list):
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)       
        x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
        cc = plt.cm.get_cmap('jet')
        cbvis=ax[qq].pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
        cax = plt.axes([0.945, 0.382, 0.01, 0.251])
        cc1=fig.colorbar(cbvis, ax=ax,cax=cax)
        cc1.set_label(label='Log$_{10}$ Viscosity (Pa $\cdot$ s)', size=23)
        cc1.ax.tick_params(labelsize=20)
        cc1.ax.yaxis.set_label_position('right')
    # ---------------------- plot setting --------------------------
        ax[qq].contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        ax[qq].set_aspect('equal')
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=26)
        ax[qq].tick_params(axis='y', labelsize=26)
        
        
        ax[qq].set_ylabel('Depth (km)',fontsize=26)
        ymajor_ticks = np.linspace(300,0,num=7)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(xmin,xmax)
        ax[qq].set_ylim(-zmin,-zmax)
        ax[-1].set_xlabel('Distance (km)',fontsize=33)
      
        # fig.savefig(figpath+model+'_vis_all.pdf')
    
if plot_pressure:
    fig, (ax)= plt.subplots(len(frame_list),1,figsize=(13,18))
    #--------------------- pressure plotting -------------------------
    for qq,frame in enumerate(frame_list):
        
        x,z,ele_x,ele_z,new_pre = dynamics_pressure(frame)
        ck = plt.cm.get_cmap('RdYlBu_r')
        # x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
        cbpre=ax[qq].pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
        cax = plt.axes([0.945, 0.385, 0.01, 0.251])
        cc1=fig.colorbar(cbpre, ax=ax,cax=cax)
        cc1.set_label(label='Pressure (MPa)', size=25)
        cc1.ax.tick_params(labelsize=20)
        cc1.ax.yaxis.set_label_position('right')
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)       
    # ---------------------- plot setting --------------------------
        ax[qq].contour(xt,-zt,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=2)
        ax[qq].set_aspect('equal')
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=26)
        ax[qq].tick_params(axis='y', labelsize=26)
        
        
        ax[qq].set_ylabel('Depth (km)',fontsize=26)
        ymajor_ticks = np.linspace(300,0,num=7)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(xmin,xmax)
        ax[qq].set_ylim(-zmin,-zmax)
        ax[-1].set_xlabel('Distance (km)',fontsize=33)
      
        # fig.savefig(figpath+model+'_pre_all.pdf')
if plot_sxx:
    fig, (ax)= plt.subplots(len(frame_list),1,figsize=(13,18))
    cn= plt.cm.get_cmap('seismic')
    #--------------------- plotting -------------------------
    for qq,frame in enumerate(frame_list):
        sxx = fl.read_sxx(frame)*100
        x,z = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)   
        cbsxx=ax[qq].pcolormesh(x,-z,sxx,cmap=cn,vmin=-200,vmax=200)
        ax[qq].set_title('sxx',fontsize=25)
        cax = plt.axes([0.13, 0.04, 0.77, 0.025])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax,orientation='horizontal')
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[qq].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
        ax[qq].set_aspect('equal')
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=23)
        ax[qq].tick_params(axis='y', labelsize=23)
        ax[qq].set_ylabel('Depth (km)',fontsize=28)
        ax[qq].set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
        ax[qq].set_ylabel('Depth (km)',fontsize=26)
        ymajor_ticks = np.linspace(300,0,num=7)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(xmin,xmax)
        ax[qq].set_ylim(-zmin,-zmax)
        ax[-1].set_xlabel('Distance (km)',fontsize=33)
        
    # fig.savefig(figpath+'Sxx_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')

if plot_srII:
    fig, (ax)= plt.subplots(len(frame_list),1,figsize=(13,18))
    cn= plt.cm.get_cmap('seismic')
    #--------------------- plotting -------------------------
    for qq,frame in enumerate(frame_list):
        sxx = fl.read_sxz(frame)*100
        # sxx = fl.read_strain(frame)[2]*100
        x,z = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)   
        cbsxx=ax[qq].pcolormesh(x,-z,sxx,cmap=cn,vmin=-100,vmax=100)
        ax[qq].set_title('sxx',fontsize=25)
        cax = plt.axes([0.13, 0.04, 0.77, 0.025])
        cc1=fig.colorbar(cbsxx, ax=ax[0],cax=cax,orientation='horizontal')
        cc1.ax.tick_params(labelsize=20)
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax[qq].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
        ax[qq].set_aspect('equal')
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=23)
        ax[qq].tick_params(axis='y', labelsize=23)
        ax[qq].set_ylabel('Depth (km)',fontsize=28)
        ax[qq].set_title('Time '+str(np.round(fl.time[frame-1],0))+' Myr',fontsize=25)
        
        ax[qq].set_ylabel('Depth (km)',fontsize=26)
        ymajor_ticks = np.linspace(300,0,num=7)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(xmin,xmax)
        ax[qq].set_ylim(-zmin,-zmax)
        ax[-1].set_xlabel('Distance (km)',fontsize=33)
        
    # fig.savefig(figpath+'Sxx_compare_of_'+str(frame)+'.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
