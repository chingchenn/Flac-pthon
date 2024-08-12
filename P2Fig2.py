#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 17:52:13 2023

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
# Model
Cocos           = 1
Nazca           = 0
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
#path = 'F:/model/'
#path = 'D:/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation/'

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

if Cocos:
    xmin,xmax = 500,900
    zmin,zmax = -150,10
    model = 'Ref_Cocos'
    model = 'Cocos_a0101'
    frame1 = 50
    frame2 = 70
    frame3 = 110
    frame4 = 180
elif Nazca:
    xmin,xmax = -100,700
    zmin,zmax = -200,10
    model = 'Nazca_aa06'
    # model = 'Nazca_v2_01'
    #model = 'Nazca_a0702'
    frame1 = 30
    frame2 = 60
    frame3 = 120
    frame4 = 180

os.chdir(path+model)
fl = flac.Flac()
time=fl.time
bwith = 3
fontsize=30
labelsize=30
### tpopgraphy
x,z = fl.read_mesh(150)
ele_x,ele_z = flac.elem_coord(x, z)
xtt = x[:,0]
ztt = z[:,0]
trench_x=xtt[np.argmin(ztt)]
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

if Cocos:
    fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(34,26))
if Nazca:
    fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(17,22))
    # fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))
amount_of_magma = 1e-6
#----------------------------- FIG1 -----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
magma_chamber = fl.read_fmagma(frame) * 100
ax.plot(ele_x[:,0]-trench_x,-ele_z[:,0],c = 'k',lw=2,linestyle='dashed')
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
if Nazca:
    file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
ax.pcolormesh(x-trench_x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax.contour(xt-trench_x,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax.scatter(ele_x[melt>1e-3]-trench_x,-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
# ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],s=100,c=melt[melt>1e-3]*1e4,cmap='gist_heat',vmin=0, vmax=100)
x_moho,z_moho=find_moho(frame,x,z,ph)

if Nazca:
    ax.text(250-trench_x,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
if Cocos:
    ax.text(520,120,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)    
#----------------------------- FIG2 -----------------------------
frame = frame2
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
magma_chamber = fl.read_fmagma(frame) * 100
ax2.plot(ele_x[:,0]-trench_x,-ele_z[:,0],c = 'k',lw=2,linestyle='dashed')
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
if Nazca:
    file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
ax2.pcolormesh(x-trench_x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax2.contour(xt-trench_x,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax2.scatter(ele_x[melt>1e-3]-trench_x,-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='#DC143C')
x_moho,z_moho=find_moho(frame,x,z,ph)
if Nazca:
    ax2.text(250-trench_x,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
if Cocos:
    ax2.text(520,120,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
#----------------------------- FIG3 -----------------------------
frame = frame3
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
magma_chamber = fl.read_fmagma(frame) * 100
ax3.plot(ele_x[:,0]-trench_x,-ele_z[:,0],c = 'k',lw=2,linestyle='dashed')
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
if Nazca:
    file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
ax3.pcolormesh(x-trench_x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax3.contour(xt-trench_x,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax3.scatter(ele_x[melt>1e-3]-trench_x,-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='#DC143C')
# ax3.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],s=100,c=melt[melt>1e-3]*1e4,cmap='gist_heat',vmin=0, vmax=100)
x_moho,z_moho=find_moho(frame,x,z,ph)

if Nazca:
    ax3.text(250-trench_x,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
if Cocos:
    ax3.text(520,120,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
#----------------------------- FIG4 -----------------------------
frame = frame4
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
melt=fl.read_fmelt(frame)
magma_chamber = fl.read_fmagma(frame) * 100
ax4.plot(ele_x[:,0]-trench_x,-ele_z[:,0],c = 'k',lw=2,linestyle='dashed')
#--------------------- phase plotting -------------------------
file=model+'_frame'+str(frame)+'_phase.grd'
if Nazca:
    file=model+'_phase_'+str(frame)+'.grd'
data = Dataset(savepath+file, mode='r')
x = data.variables['x'][:]
z = data.variables['y'][:]
ph = data.variables['z'][:]
ax4.pcolormesh(x-trench_x,-z,ph,cmap=phase8,vmin=1, vmax=20)
ax4.contour(xt-trench_x,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
ax4.scatter(ele_x[melt>1e-3]-trench_x,-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='#DC143C')
# ax4.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],s=100,c=melt[melt>1e-3]*1e4,cmap='gist_heat',vmin=0, vmax=100)
x_moho,z_moho=find_moho(frame,x,z,ph)

if Nazca:
    ax4.text(250-trench_x,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
if Cocos:
    ax4.text(520,120,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
# ---------------------- plot setting --------------------------
if Nazca:
    xmajor_ticks=np.array([-100,0,100,200,300,400,500,600,700])
if Cocos:
    xmajor_ticks = np.linspace(500,900,num=5)
ymajor_ticks = np.linspace(200,0,num=5)
for aa in [ax,ax2,ax3,ax4]:#,ax5,ax6,ax7,ax8]:
    aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
    aa.set_aspect('equal')
    aa.set_yticks(ymajor_ticks)
    aa.set_ylabel('depth (km)',fontsize=fontsize)
    aa.set_xticks(xmajor_ticks)
    aa.set_xlim(xmin,xmax)
    aa.set_ylim(-zmin,-zmax)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
ax4.set_xlabel('distance (km)',fontsize=fontsize)
if png:
    fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase_magmatism.png')
if pdf:
    fig.savefig(figpath+'figure1_v7.pdf')