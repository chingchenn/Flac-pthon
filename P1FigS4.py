#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 13:06:54 2023

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


plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- DO WHAT -----------------------------------
### pdf or png
png             = 0
pdf             = 0

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'

phase_uppercrust = 2
phase_basalt = 3
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

xmin,xmax = 250,1000
zmin,zmax = -200,10
model = 'Nazca_aa06'
#model = 'Nazca_a0702'
frame1 = 30
frame2 = 60
frame3 = 140
frame4 = 170

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

def find_moho(frame,x,z,ph):
    x_moho=np.zeros(len(x))
    z_moho=np.zeros(len(x))
    for kk in range(len(x)-1,0,-1):
        bb = z[np.where(ph[:,kk]==phase_basalt)]
        ee = z[np.where(ph[:,kk]==phase_eclogite)]
        if len(bb)==0 and len(ee) ==0:
            continue
        x_moho[kk]=x[kk]
        if len(bb)==0:
            z_moho[kk] = np.min(ee)
        elif len(ee)==0:
            z_moho[kk] = np.min(bb)
        else:
           z_moho[kk]=min(np.min(bb),np.min(ee))
    return x_moho, z_moho
#------------------------------------------------------------------------------
fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,22))
#----------------------------- FIG1 -----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
phh=ph.data[ph.data>0]
ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
ax.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
#----------------------------- FIG2 -----------------------------
frame = frame2
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
phh=ph.data[ph.data>0]
ax2.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax2.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
ax2.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
#----------------------------- FIG3 -----------------------------
frame = frame3
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
phh=ph.data[ph.data>0]
ax3.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax3.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
ax3.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)

#----------------------------- FIG4 -----------------------------
frame = frame4
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
phh=ph.data[ph.data>0]
ax4.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax4.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
ax4.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)

# ---------------------- plot setting --------------------------
xmajor_ticks=np.array([300,400,500,600,700,800,900,1000])
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
    fig.savefig(figpath+'FigS4v2.pdf')