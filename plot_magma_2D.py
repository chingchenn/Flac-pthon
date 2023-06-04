#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 22:50:33 2022

@author: chingchen
"""
import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
from numpy import unravel_index
import matplotlib.colors as mcolors

path = '/home/jiching/geoflac/'
#path = '/home/jiching/test_geoflac/geoflac/'
#path = '/home/jiching/geoflac_T/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/04model/'
#path = '/scratch2/jiching/'
path='/Users/chingchen/Desktop/model/'
#savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
savepath='/Users/chingchen/Desktop/data/'
#figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
figpath='/Users/chingchen/Desktop/figure/'
# model = sys.argv[1]
fig1=1 # Nazca 
fig2=0 # Cocos
fig3=0 # Cocos without basalt-to-eclogite


plt.rcParams["font.family"] = "Times New Roman"
bwith = 3

phase_uppercrust = 2
phase_oceanic = 3
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

model='Nazca_a0702'
os.chdir(path+model)
fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
# end = 250
time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
rainbow = cm.get_cmap('gray_r',end)
meltcolor = cm.get_cmap('turbo',end)
newcolors = rainbow(np.linspace(0, 1, end))
time_color = meltcolor(np.linspace(0,1,end))
if fig1:
    fig, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(15,12),gridspec_kw={'height_ratios':[3,3,2]})
    xxx_trench = np.max(x_trench)
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) * 100
        melt = fl.read_fmelt(i) * 100
        ax2.scatter(ele_x[magma_chamber>1.5e-3],-ele_z[magma_chamber>1.5e-3],color=time_color[i],zorder=1,s=10)
        # time = fl.time[i]
        qqq=ax1.scatter(ele_x[melt>0.1],-ele_z[melt>0.1],color=time_color[i],s = 10)
    def x2dis(x):
        return x -xxx_trench
    def dis2x(x):
        return x +xxx_trench
    ax3 = ax1.secondary_xaxis('top', functions=(x2dis,dis2x))
    ax4 = ax2.secondary_xaxis('top', functions=(x2dis,dis2x))
    #ax3.set_xlabel('Distance from trench (km)',fontsize=30)
    
    for ax in [ax1,ax2,ax5]:
        ax.grid()
        ax.set_xlim(300,1000)
        ax.set_ylim(150,0)
        ax.set_ylabel('Depth (km)',fontsize=26)
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.set_aspect('equal')
    for ax in [ax3,ax4]:
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
    
    for frame in [10.0, 20.0, 30.0, 40.0]:
        xx,zz = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
        xx=xx[zz<0]
        zz=zz[zz<0]
        ax1.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        
        
    x_change = np.zeros(end)
    z_change = np.zeros(end)
    time = np.ones(end)
    for i in range(2,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        phase = fl.read_phase(i)
        time[i] = i * 0.2
        for xx in range(len(ele_x)):
            if True in (phase[xx,:]==phase_eclogite):
                for zz in range(len(ele_x[0])):
                    if (phase[xx,zz]==phase_eclogite):
                        x_change[i]=ele_x[xx,zz]
                        z_change[i]=-ele_z[xx,zz]
                        break
                break
    cbtime=ax5.scatter(x_change[x_change>0],z_change[x_change>0],c=time[x_change>0],cmap = 'turbo',vmin=0,vmax=40)
    ax5.set_ylim(100,0)
    ax5.set_title('Location of basalt to eclogite phase change',fontsize=28)
    ax5.set_xlabel('Distance (km)',fontsize=28)
    ax1.set_title('Location of partial melting',fontsize=28)
    ax2.set_title('Location of magma chamber',fontsize=28)
    cax = plt.axes([0.99, 0.082, 0.03, 0.850])
    cc1=fig.colorbar(cbtime, ax=ax,cax=cax)
    cc1.set_label(label='Time (Myr)', size=23)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('right')
    fig.tight_layout()
    # fig.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Ref_Nazca/'+'Nazca_a0702_2Dtime_series_7.pdf')



model='Ref_Cocos'
os.chdir(path+model)
fl = flac.Flac()
end = fl.nrec
end = 200
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
rainbow = cm.get_cmap('gray_r',end)
meltcolor = cm.get_cmap('rainbow',end)
newcolors = rainbow(np.linspace(0, 1, end))
time_color = meltcolor(np.linspace(0,1,end))
if fig2:
    fig2, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(13,14),gridspec_kw={'height_ratios':[3,3,2]})
    xxx_trench = np.average(x_trench)
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) * 100
        melt = fl.read_fmelt(i) * 100
        ax2.scatter(ele_x[magma_chamber>1e-1],-ele_z[magma_chamber>1e-1],color=time_color[i],zorder=1,s=10)
        # time = fl.time[i]
        qqq=ax1.scatter(ele_x[melt>1e-2],-ele_z[melt>1e-2],color=time_color[i],s = 10)
    def x2dis(x):
        return x -xxx_trench
    def dis2x(x):
        return x +xxx_trench
    ax3 = ax1.secondary_xaxis('top', functions=(x2dis,dis2x))
    ax4 = ax2.secondary_xaxis('top', functions=(x2dis,dis2x))
    # ax3.set_xlabel('Distance from trench (km)',fontsize=30)
    
    for ax in [ax1,ax2,ax5]:
        ax.grid()
        ax.set_xlim(500,900)
        ax.set_ylim(150,0)
        ax.set_ylabel('Depth (km)',fontsize=26)
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.set_aspect('equal')
    for ax in [ax3,ax4]:
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
    
    for frame in [10.0, 20.0, 30.0, 40.0]:
        xx,zz = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
        xx=xx[zz<0]
        zz=zz[zz<0]
        ax1.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
    
    x_change = np.zeros(end)
    z_change = np.zeros(end)
    time = np.ones(end)
    for i in range(2,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        phase = fl.read_phase(i)
        time[i] = i * 0.2
        for zz in range(len(ele_x[0])):
            if True in (phase[:,zz]==phase_eclogite):
                for xx in range(len(ele_x)):
                    if ele_x[xx,zz]<700:
                        # print(z_change[i]- z_change[i-1])
                        if (phase[xx,zz]==phase_eclogite) and (-ele_z[xx,zz]> z_change[i-1]):
                            x_change[i]=ele_x[xx,zz]
                            z_change[i]=-ele_z[xx,zz]
                            break
                    else:
                        if (phase[xx,zz]==phase_eclogite):
                            x_change[i]=ele_x[xx,zz]
                            z_change[i]=-ele_z[xx,zz]
                            break
                break
    cbtime=ax5.scatter(x_change[x_change>0],z_change[x_change>0],c=time[x_change>0],cmap = 'rainbow',vmin=0,vmax=40)
    ax5.set_ylim(100,0)
    ax5.set_title('Location of basalt to eclogite phase change',fontsize=28)
    ax1.set_title('Location of partial melting',fontsize=28)
    ax2.set_title('Location of magma chamber',fontsize=28)
    ax5.set_xlabel('Distance (km)',fontsize=28)
    cax = plt.axes([0.945, 0.137, 0.03, 0.725])
    cc1=fig2.colorbar(cbtime, ax=ax,cax=cax)
    cc1.set_label(label='Time (Myr)', size=23)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    fig2.tight_layout()
    # fig2.savefig(figpath+model+'_2Dtime_series.pdf')
if fig3:
    fig3, (ax1,ax2)= plt.subplots(2,1,figsize=(13,12))
    xxx_trench = np.average(x_trench)
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) * 100
        melt = fl.read_fmelt(i) * 100
        ax2.scatter(ele_x[magma_chamber>1e-1],-ele_z[magma_chamber>1e-1],color=time_color[i],zorder=1,s=10)
        # time = fl.time[i]
        qqq=ax1.scatter(ele_x[melt>1e-2],-ele_z[melt>1e-2],color=time_color[i],s = 10)
    def x2dis(x):
        return x -xxx_trench
    def dis2x(x):
        return x +xxx_trench
    ax3 = ax1.secondary_xaxis('top', functions=(x2dis,dis2x))
    ax4 = ax2.secondary_xaxis('top', functions=(x2dis,dis2x))
    # ax3.set_xlabel('Distance from trench (km)',fontsize=30)
    
    for ax in [ax1,ax2]:
        ax.grid()
        ax.set_xlim(500,900)
        ax.set_ylim(150,0)
        ax.set_ylabel('Depth (km)',fontsize=26)
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.set_aspect('equal')
    for ax in [ax3,ax4]:
        ax.tick_params(axis='x', labelsize=26 )
        ax.tick_params(axis='y', labelsize=26 )
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
    
    for frame in [10.0, 20.0, 30.0, 40.0]:
        xx,zz = np.loadtxt(savepath+model+'_'+str(frame)+'_final_slab.txt').T
        xx=xx[zz<0]
        zz=zz[zz<0]
        xx = fd.moving_window_smooth(xx, 3)
        zz = fd.moving_window_smooth(zz, 3)
        qqq=ax1.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
    
    ax2.set_xlabel('Distance (km)',fontsize=28)
    ax1.set_title('Location of partial melting',fontsize=28)
    ax2.set_title('Location of magma chamber',fontsize=28)
    cax = plt.axes([1.03, 0.137, 0.03, 0.7])
    norm1 = mcolors.Normalize(vmin=0, vmax=40)
    cc1=fig3.colorbar(cm.ScalarMappable(norm=norm1,cmap='rainbow'), ax=ax,cax=cax)
    cc1.set_label(label='Time (Myr)', size=23)
    cc1.ax.tick_params(labelsize=20)
    cc1.ax.yaxis.set_label_position('left')
    fig3.tight_layout()
    # fig2.savefig(figpath+model+'_2Dtime_series.pdf')
