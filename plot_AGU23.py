#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:35:47 2023

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
Cocos           = 0
Nazca           = 1
### pdf or png
png             = 0
pdf             = 0
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/AGU2023FLAC/'

plot_rock_field = 1
nazca_magma = 0 ## Nazca melting phase and melting location v.s. trench
cocos_magma = 0 ## Cocos melting phase and melting location v.s. trench
cocos_magma_molten = 0
cocos_magma_molten_save = 0
magma_inone = 0
nazca_magma_save = 0
cocos_magma_save = 0
magma_save  = 0
zoomin_melting_mex = 0
zoomin_melting_mex_save = 0
plot_2dmagma_nazca = 0
plot_2dmagma_nazca_save = 0
plot_s1_profile_nazca = 0
plot_s1_profile_cocos = 0
plot_sII_profile_nazca = 0
plot_sII_profile_cocos = 0



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
    frame1 = 50
    frame2 = 70
    frame3 = 110
    frame4 = 180
elif Nazca:
    xmin,xmax = 250,1000
    zmin,zmax = -200,10
    model = 'Nazca_aa06'
    #model = 'Nazca_a0702'
    frame1 = 30
    frame2 = 60
    frame3 = 120
    frame4 = 178

os.chdir(path+model)
fl = flac.Flac()
time=fl.time
bwith = 2
fontsize=30
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
def get_magma(end):
    #-----------Time Series-----------
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    rrr=np.zeros(end)
    depth_magma = np.zeros(end)
    for i in range(1,end):
        x,z = fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        chamber=fl.read_fmagma(i)
        melt[i] = np.max(mm)
        magma[i] = np.max(chamber)
        ele_x,ele_z = flac.elem_coord(x, z)
        ui,uj= np.unravel_index(chamber.argmax(), chamber.shape)
        depth_magma[i] = -ele_z[ui,uj]
        if  magma[i]!=0:
            rrr[i]= melt[i]/magma[i]
        arc_vol[i]=np.sum(fl.read_area(i)[phase ==14])/1e6
        yymelt[i]=(fl.read_fmelt(i)*fl.read_area(i)/1e6).sum()
        yychamber[i]=(fl.read_fmagma(i)*fl.read_area(i)/1e6).sum()
    return melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma
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

colors2=[
 '#C98F49', '#92C0DF', '#2553C7', '#FFFFFF', '#6495ED',
 '#2E8B57', '#524B52', '#9A32CD', '#6B8E23','#D4DBF4',
 '#D8BFD8','#999999','#F2C85B','#92C0DF','#999999',
 '#4CC552','#999999','#999999','#999999','#999999']
phase8= matplotlib.colors.ListedColormap(colors2)

if plot_rock_field:
    if Cocos:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(14,20))
    if Nazca:
        fig, (ax,ax2,ax3,ax4)= plt.subplots(4,1,figsize=(34,22))
    amount_of_magma = 1e-6
    #----------------------------- FIG1 -----------------------------
    frame = frame1
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    magma_chamber = fl.read_fmagma(frame) * 100
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    if Nazca:
        file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4,c='#DC143C')
    x_moho,z_moho=find_moho(frame,x,z,ph)  
    #----------------------------- FIG2 -----------------------------
    
    frame = frame2
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    magma_chamber = fl.read_fmagma(frame) * 100
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    if Nazca:
        file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax2.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax2.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    x_moho,z_moho=find_moho(frame,x,z,ph)
    #----------------------------- FIG3 -----------------------------
    frame = frame3
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    magma_chamber = fl.read_fmagma(frame) * 100
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    if Nazca:
        file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax3.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax3.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    x_moho,z_moho=find_moho(frame,x,z,ph)
    #----------------------------- FIG4 -----------------------------
    frame = frame4
    xt,zt = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(xt, zt)
    temp = fl.read_temperature(frame)
    melt=fl.read_fmelt(frame)
    magma_chamber = fl.read_fmagma(frame) * 100
    #--------------------- phase plotting -------------------------
    file=model+'_frame'+str(frame)+'_phase.grd'
    if Nazca:
        file=model+'_phase_'+str(frame)+'.grd'
    data = Dataset(savepath+file, mode='r')
    x = data.variables['x'][:]
    z = data.variables['y'][:]
    ph = data.variables['z'][:]
    ax4.pcolormesh(x,-z,ph,cmap=phase8,vmin=1, vmax=20)
    ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
    ax4.scatter(ele_x[melt>1e-3],-ele_z[melt>1e-3],melt[melt>1e-3]*1e4*3,c='#DC143C')
    x_moho,z_moho=find_moho(frame,x,z,ph)
    # ---------------------- plot setting --------------------------
    if Nazca:
        xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])
    if Cocos:
        xmajor_ticks = np.linspace(500,900,num=5)
    ymajor_ticks = np.linspace(200,0,num=5)
    for aa in [ax,ax2,ax3,ax4]:
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
        fig.savefig(figpath+'Mexico_phase_field.png')
    if pdf:
        fig.savefig(figpath+'Mexico_phase_field.pdf')


if nazca_magma:
    fig6, (ax2,ax1)= plt.subplots(2,1,figsize=(15,12),gridspec_kw={'height_ratios':[1,0.5]})  
    model='Nazca_aa06'
    model='Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    ####====================== FIG1 melting phase ======================
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    ####====================== FIG2 melting location ======================
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=70,cmap='OrRd',vmax=1,vmin=-0)
        hod = ele_x[melting>6e-3]
        hhh = (melting>6e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=70,cmap='OrRd',vmax=1,vmin=-0)
    #melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma=get_magma(end)
    cax = plt.axes([0.999, 0.42, 0.02, 0.56])
    cbar=fig6.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('melting percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    #================================figure setting================================
    name='Nazca_a0702'+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_ylim(0,10)    
    ax2.set_ylim(0,600)
    ax2.set_ylabel('distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
    ax1.set_xlabel('time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    if nazca_magma_save:
        fig6.savefig(figpath+'melting_nazca.pdf')
if cocos_magma:
    fig6, (ax2,ax1)= plt.subplots(2,1,figsize=(18,12),gridspec_kw={'height_ratios':[1,0.5]})  
    import matplotlib.gridspec as gridspec

    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(gs[0, 0]) # row 0, col 0
    ax1.plot([0,1])
    
    ax2 = fig.add_subplot(gs[0, 1]) # row 0, col 1
    ax2.plot([0,1])
    
    ax3 = fig.add_subplot(gs[1, :]) # row 1, span all columns
    ax3.plot([0,1])
    model='Ref_Cocos'
    #model='Cocos_aa02'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    ####====================== FIG1 melting phase ======================
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    ####====================== FIG2 melting location ======================
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=70,cmap='OrRd',vmax=1,vmin=-0)
        hod = ele_x[melting>6e-3]
        hhh = (melting>6e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=70,cmap='OrRd',vmax=1,vmin=-0)

    cax = plt.axes([0.999, 0.42, 0.02, 0.56])
    cbar=fig6.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('melting percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')

    #================================figure setting================================
    name=model+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        #aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_ylim(0,4)    
    ax2.set_ylim(0,300)
    ax2.set_ylabel('distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
    ax1.set_xlabel('time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    # ax1.legend(fontsize=fontsize,facecolor='white')
    if cocos_magma_save:
        fig6.savefig(figpath+'melting_Mexico_time.pdf')
if cocos_magma_molten:
    fig7, (ax1,ax2)= plt.subplots(1,2,figsize=(22,5),gridspec_kw={'width_ratios':[15,27]})  
    ####====================== FIG1 melting phase ======================
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    ####====================== FIG2 melting phase ======================
    axx=ax2.twinx()
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    kkk2 = fd.moving_window_smooth(phase_p4[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    # ax3.bar(time[total>0],yychamber[total>0],width=0.17,color='orange',alpha=0.5)
    axx.plot(time[total>0],kkk,c='#4169E1')
    axx.plot(time[total>0],kkk2,c='seagreen')
    ax2.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax2.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax2.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    #================================figure setting================================
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.set_xlabel('time (Myr)',fontsize=fontsize)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_xlim(0,15)    
    ax2.set_xlim(13,40)
    ax2.set_ylim(0,1)
    axx.tick_params(labelsize=fontsize,labelcolor="#8A2BE2")
    axx.set_ylim(-30,70)
    axx.set_ylabel('% of melting rocks',fontsize=fontsize,color ="#8A2BE2" )
    ax1.set_ylabel('meling rocks (km$^3$/km)',fontsize=fontsize)
    if cocos_magma_molten_save:
        fig7.savefig(figpath+'molten_Mexico_time.pdf')

if magma_inone:
    fig6, (ax1,ax2,ax3)= plt.subplots(3,1,figsize=(20,18),gridspec_kw={'height_ratios':[1,1,1]})  
    model='Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    #----------------------------------- FIG1 melting phase -----------------------------------
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    name=model+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    ax1.vlines(x=time_flat[0]-0.5, ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
    ##----------------------------------- FIG3 melting location -----------------------------------
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        ttt = np.ones(len(hod))*ii*0.2
        ax3.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c ='#708090',s=70)
    ax3.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c ='#708090',s=70, label = 'chilean')
    #----------------------------------- FIG2 melting phase -----------------------------------
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax2.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax2.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax2.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    #---------------------------------- FIG3 melting location -----------------------------------
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(3,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        ttt = np.ones(len(hod))*ii*0.2
        ax3.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c ='#AE6378' ,s=70)
    ax3.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c ='#AE6378' ,s=70, label = 'mexican')
    name=model+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    ax2.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
    #---------------------------------- figure setting -----------------------------------
    for aaa in [ax2,ax1,ax3]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_ylim(0,10) 
    ax2.set_ylim(0,4)    
    ax3.set_ylim(0,600)
    ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
    ax2.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
    ax3.set_ylabel('distance from trench (km)',fontsize=fontsize)
    ax3.set_xlabel('time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    ax3.legend(fontsize=fontsize,facecolor='white')
    if magma_save:
        fig6.savefig(figpath+'melting.pdf')
if zoomin_melting_mex:
    fig7, (ax1)= plt.subplots(1,1,figsize=(18,5))  
    model='Ref_Cocos'
    os.chdir(path+model)
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    #================================figure setting================================
    for aaa in [ax1]:
        aaa.tick_params(axis='x', labelsize=fontsize-5)
        aaa.tick_params(axis='y', labelsize=fontsize-5)
        aaa.set_xlim(13.4,40)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith-1)
    ax1.set_ylim(0,0.6)
    # ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    # ax1.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig7.tight_layout()
    # ax1.legend(fontsize=fontsize,facecolor='white')
    if zoomin_melting_mex_save:
        fig7.savefig(figpath+'melting_Mexico_time_enlarge.pdf')
   # 
if plot_2dmagma_nazca:
    # phase_uppercrust = 2
    # phase_oceanic = 3
    # phase_mantle1 = 4
    # phase_schist = 5
    # phase_mantle2 = 8
    # phase_serpentinite = 9
    # phase_sediment = 10
    # phase_sediment_1 = 11
    # phase_eclogite = 13
    # phase_lowercrust = 14
    # phase_hydratedmantle = 16
    # phase_oceanic_1 = 17
    # phase_eclogite_1 = 18

    model='Nazca_aa06'
    os.chdir(path+model)
    fl = flac.Flac()
    end = fl.nrec
    nex = fl.nx - 1
    nez = fl.nz - 1
    time = fl.time
    cmap = 'jet'

    time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
    rainbow = cm.get_cmap('gray_r',end)
    meltcolor = cm.get_cmap('turbo',end)
    newcolors = rainbow(np.linspace(0, 1, end))
    time_color = meltcolor(np.linspace(0,1,end))

    fig8, (ax2,ax5)= plt.subplots(2,1,figsize=(20,8),gridspec_kw={'height_ratios':[2,0.15]})
    xxx_trench = np.max(x_trench)
    chamber_limit = 1e-3
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) 
        melt = fl.read_fmelt(i) * 100
        ax2.scatter(ele_x[magma_chamber>chamber_limit],-ele_z[magma_chamber>chamber_limit],c=time_color[i],cmap =cmap,s=10,vmin=0,vmax=40)
    for ax in [ax2]:
        ax.grid()
        ax.set_xlim(300,1000)
        ax.set_ylim(150,0)
        ax.set_ylabel('depth (km)',fontsize=fontsize)
        ax.tick_params(labelsize=fontsize )
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(bwith)
        ax.set_aspect('equal')
        ymajor_ticks = np.linspace(150,0,num=4)
        ax.set_yticks(ymajor_ticks)
    for frame in [10.0, 20.0, 30.0, 40.0]:
        xx,zz,xt = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
        xx=xx[zz<0]
        zz=zz[zz<0]
        ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        
    ax2.set_xlabel('distance (km)',fontsize=fontsize)
    # ax1.set_title('location of partial melting',fontsize=fontsize)
    # ax2.set_title('location of magma fraction',fontsize=fontsize)
    norm = mpl.colors.Normalize(vmin=0,vmax=40)
    cc1 = mpl.colorbar.ColorbarBase(ax5,cmap=cmap,norm = norm, orientation='horizontal')
    cc1.set_label(label='time (Myr)', size=fontsize)
    cc1.ax.tick_params(labelsize=fontsize)
    if plot_2dmagma_nazca_save:
        fig8.savefig(figpath+'Nazca_magmma2D.pdf')


if plot_s1_profile_nazca:  
    model='Nazca_aa06'
    shift = 320
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    frame = 200
    skip = (slice(None, None, 5), slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    fontsize=33
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('seismic_r')
    cbsxx = plt.cm.get_cmap('bwr_r')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx*100,cmap=cbsxx,vmin=-500,vmax=500,shading='gouraud')
    
    ax.quiver(ele_x[skip],-ele_z[skip],s1.T[skip],s3.T[skip],pivot = 'mid',headwidth=0,color = 'green')#,
              #angles='xy', scale_units='xy', scale=0.1,)
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
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(230,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(300,950)
    
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
    
    ax.set_xlabel('distance (km)',fontsize=fontsize)
    ax.set_ylabel('depth (km)',fontsize=fontsize)
    # fig.savefig(figpath+model+'nazca_sxx_gourand.pdf')

if plot_s1_profile_cocos:
    model = 'Ref_Cocos'
    shift = 550
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    frame = 200
    skip = (slice(None, None, 5), slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    fontsize=36
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('seismic_r')
    cbsxx = plt.cm.get_cmap('bwr_r')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx*100,cmap=cbsxx,vmin=-500,vmax=500,shading='gouraud')
    
    ax.quiver(ele_x[skip],-ele_z[skip],s1.T[skip],s3.T[skip],pivot = 'mid',headwidth=0,color = 'green',)
              # angles='xy', scale_units='xy', scale=0.1,)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=fontsize)
    ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(150,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(500,850)
    #xmajor_ticks = np.linspace(250,1000,num=7)
    #ax.set_xticks(xmajor_ticks)
    
   
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Ref_Cocos_30.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=5)
       
    xxm,zzm,xtm = np.loadtxt(savepath+'Ref_Cocos_30_final_moho_slab.txt').T
    xxmm=xxm[(zzm<0)*(xxm>0)]
    zzmm=zzm[(zzm<0)*(xxm>0)]
    zzm = fd.moving_window_smooth(zzm,3)
    ax.plot(xxmm+shift,-zzmm+2,color='k',lw=5)
    
    
    ax.set_xlabel('distance (km)',fontsize=fontsize)
    ax.set_ylabel('depth (km)',fontsize=fontsize)
    # fig.savefig(figpath+model+'cocos_sxx_gourand.pdf')
    

if plot_sII_profile_nazca:
    model='Nazca_aa06'
    shift = 320
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    frame = 200
    skip = (slice(None, None, 5), slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    fontsize=33
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    sII = fl.read_sII(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('hot_r')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sII*100,cmap=cbsxx,vmin=0,vmax=1000,shading='gouraud')
    magnitude = np.sqrt(s1.T[skip]**2+s3.T[skip]**2)
    vector = ax.quiver(ele_x[skip][(magnitude>1)],-ele_z[skip][(magnitude>1)],
              s1.T[skip][(magnitude>1)],s3.T[skip][(magnitude>1)],
              pivot = 'mid',headwidth=0,color = 'green',
                scale_units='xy', scale=0.08,)
    ax.quiverkey(vector,350,170,5,"5 cm/y",coordinates='data',color='k',fontproperties={'size': labelsize-3})
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    cc1.set_label(label='sII (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=1.5)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(230,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(300,950)
    
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=2)
       
    xxm,zzm,xtm = np.loadtxt(savepath+'Nazca_aa06_40_final_moho_slab.txt').T
    xxm=xxm[zzm<0]
    zzm=zzm[zzm<0]
    zzm = fd.moving_window_smooth(zzm,3)
    ax.plot(xxm+shift-5,-zzm+2,color='k',lw=2)
    
    ax.set_xlabel('distance (km)',fontsize=fontsize)
    ax.set_ylabel('depth (km)',fontsize=fontsize)
    fig.savefig(figpath+model+'nazca_sII_gourand.pdf')
if plot_sII_profile_cocos:
    model = 'Ref_Cocos'
    shift = 550
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    frame = 200
    skip = (slice(None, None, 5),slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    fontsize=36
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    sII = fl.read_sII(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('hot_r')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sII*100,cmap=cbsxx,vmin=0,vmax=1000,shading='gouraud')
    magnitude = np.sqrt(s1.T[skip]**2+s3.T[skip]**2)
    vector = ax.quiver(ele_x[skip][(magnitude>1)],-ele_z[skip][(magnitude>1)],
              s1.T[skip][(magnitude>1)],s3.T[skip][(magnitude>1)],
              pivot = 'mid',headwidth=0,color = 'green',
                scale_units='xy', scale=0.08,)
    ax.quiverkey(vector,550,130,5,"5 cm/y",coordinates='data',color='k',fontproperties={'size': labelsize-3})
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=fontsize)
    cc1.set_label(label='sII (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=1.5)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(150,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(500,850)
    #xmajor_ticks = np.linspace(250,1000,num=7)
    #ax.set_xticks(xmajor_ticks)
    
   
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Ref_Cocos_30.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=2)
       
    xxm,zzm,xtm = np.loadtxt(savepath+'Ref_Cocos_30_final_moho_slab.txt').T
    xxmm=xxm[(zzm<0)*(xxm>0)]
    zzmm=zzm[(zzm<0)*(xxm>0)]
    zzm = fd.moving_window_smooth(zzm,3)
    ax.plot(xxmm+shift,-zzmm+2,color='k',lw=2)
    
    
    ax.set_xlabel('distance (km)',fontsize=fontsize)
    ax.set_ylabel('depth (km)',fontsize=fontsize)
    fig.savefig(figpath+model+'cocos_sII_gourand.pdf')
    