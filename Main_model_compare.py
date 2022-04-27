#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 13:11:24 2021
@author: jiching
"""
import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt

#---------------------------------- DO WHAT -----------------------------------
trench_plot             = 0
dip_plot                = 0
plate_geometry          = 1
force_plot_LR           = 0
force_plot_RF           = 0
vel_plot                = 0
stack_topo_plot         = 0
flat_slab_plot          = 0
magma_plot   	    	= 0

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
path = '/scratch2/jiching/03model/'
#path = 'F:/model/'

savepath='/home/jiching/geoflac/data/'
savepath = '/Users/ji-chingchen/Desktop/data/'
savepath = 'D:\\OneDrive - 國立台灣大學/resarch/data/'
#savepath='D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
model_list=['h0409','h0408','h0405','h0406']
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
plt.rcParams["font.family"] = "Times New Roman"
##------------------------------------ plot -----------------------------------
if trench_plot:
    print('--- start plotting the trench and topography with time ---')
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    for kk,model in enumerate(model_list):
        name='trench_for_'+model
        df = pd.read_csv(savepath+name+'.csv')
        dis,time,topo=f2.get_topo()
        ax.plot(df.trench_x[df.trench_x>0],df.time[df.trench_x>0],lw=2,label=model,color=newcolors[kk])
    ax.set_xlim(0,dis[-1][-1])
    ax.set_ylim(0,df.time[-1])
    ax.set_title(str(model)+" Bathymetry Evolution",fontsize=24)
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    fig.savefig(figpath+'compare_trench_location_'+model_list[0]+'_'+model_list[-1]+'.png')
    print('=========== DONE =============')
if dip_plot:
    print('--- start plot dip with time ---')
    fig, (ax2)= plt.subplots(1,1,figsize=(10,7))
    for kk,model in enumerate(model_list):
        name = 'plate_dip_of_'+model
        df = pd.read_csv(savepath+name+'.csv')
        ax2.plot(df.time[df.angle>0],df.angle[df.angle>0],lw=2,label=model,color=newcolors[kk])
    #ax2.set_xlim(0,df.time[-1])
    ax2.set_title('Angle Variation',fontsize=24)
    ax2.set_xlabel('Time (Myr)',fontsize=20)
    ax2.set_ylabel('Angel ($^\circ$) ',fontsize=20)
    ax2.grid()
    ax2.legend(fontsize=16)
    bwith = 3
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+'compare_dip_'+model_list[0]+'_'+model_list[-1]+'.png')
    print('=========== DONE =============')
if plate_geometry:
    print('--- start plot geometry ---')
    fig2, (ax2) = plt.subplots(1,1,figsize=(14,8))
    for kk,model in enumerate(model_list):
        xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
        xx= fd.moving_window_smooth(xmean[xmean>0], 12)
        ztop = fd.moving_window_smooth(ztop[xmean>0], 12)
        xmean=xx
        ax2.plot(xmean,ztop,c=newcolors[kk],label=model,lw=5)
    #ax2.set_xlim(0,max(xmean)+10)
    # ax2.set_title("slab comparation",fontsize=16)
    # ax2.set_ylabel("Depth (km)",fontsize=16)
    # ax2.set_xlabel("Distance relative to trench (km)",fontsize=16)
    # ax2.legend(fontsize=16)
    bwith = 3
    ax2.set_ylim(-150,0)
    ax2.set_xlim(0,500)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.set_aspect('equal')
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    fig2.savefig('D:\\OneDrive - 國立台灣大學/master03/Seminar/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.pdf')
    fig2.savefig(figpath+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.png')
    print('=========== DONE =============')
if force_plot_LR:
    print('--- start plot left and right force with time ---')
    fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
    for kk,model in enumerate(model_list):
        filepath = savepath+model+'_forc.txt'
        temp1=np.loadtxt(filepath)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        ax.scatter(time,forc_l+forc_r,s=4,label=model,color=newcolors[kk])
        ax2.scatter(time,forc_r,s=4,label=model,color=newcolors[kk])
    ax.set_xlim(0,time[-1])
    ax.set_title('oceanic side force',fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid()
    ax2.set_xlim(0,time[-1])
    ax2.set_title('continental side force',fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.legend(fontsize=16)
    fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_forc.png')
    print('=========== DONE =============')
if force_plot_RF:
    print('--- start plot ringforce with time ---')
    fig2, (ax3)= plt.subplots(1,1,figsize=(10,8))   
    for kk,model in enumerate(model_list):
        filepath =savepath+model+'_forc.txt' 
        temp1=np.loadtxt(filepath)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        ax3.scatter(time,ringforce,s=2,label=model,color=newcolors[kk])
    ax3.set_xlim(0,time[-1])
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax3.grid()
    bwith = 3
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    fig2.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_ringforc.png')
    print('=========== DONE =============')
if vel_plot:
    print('--- start plot velocity with time ---')
    fig3, (ax4)= plt.subplots(1,1,figsize=(10,8))
    for kk,model in enumerate(model_list):
        filepath = savepath+model+'_forc.txt'
        temp1=np.loadtxt(filepath)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        ax4.plot(time,vl*31545741325,lw=2,label=model,color=newcolors[kk])
    ax4.set_xlim(0,time[-1])
    ax4.set_title('oceanic side velocity',fontsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.grid()
    ax4.set_xlabel('Time (Myr)',fontsize=16)
    ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
    bwith = 3
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    fig3.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_vel.png')
    print('=========== DONE =============')
if stack_topo_plot:
    print('--- start plot topo with time ---')
    fig2, (ax2) = plt.subplots(1,1,figsize=(14,10))
    for kk,model in enumerate(model_list):
        name=model+'_stack_topography.txt'
        xmean,ztop=np.loadtxt(savepath+name).T
        ax2.plot(xmean,ztop,lw=5,label=model,color=newcolors[kk])
    ax2.set_xlabel('Distance (km)',fontsize=16)
    ax2.set_ylabel('Height (km)',fontsize=16)
    ax2.set_title('Topography',fontsize=16)
    ax2.legend(fontsize=16)
    bwith = 3
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    # fig2.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_topo_analysis.png')
    print('=========== DONE =============')
if flat_slab_plot:
    print('-----plotting flatslab-----')
    fig2, (ax1,ax2) = plt.subplots(2,1,figsize=(10,6))
    for kk,model in enumerate(model_list):
       name=model+'_flatslab_time_len.txt'
       time,length,depth=np.loadtxt(savepath+name).T
       ax1.plot(time,length,lw=3,label=model,color=newcolors[kk])
       ax2.plot(time,depth,lw=3,label=model,color=newcolors[kk])
    ax1.set_xlim(0,time[-1])
    ax1.set_title('flat slab properties',fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.grid()
    ax1.set_ylabel('length (km)',fontsize=16)
    ax2.set_xlim(0,time[-1])
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    ax2.set_xlabel('Time (Myr)',fontsize=16)
    ax2.set_ylabel('depth (km) ',fontsize=16)
    ax2.legend(fontsize=16)
    bwith = 3
    ax1.spines['bottom'].set_linewidth(bwith)
    ax1.spines['top'].set_linewidth(bwith)
    ax1.spines['right'].set_linewidth(bwith)
    ax1.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig2.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_flatslab_length.png')
    print('=========== DONE =============')
if magma_plot:
    print('--------plotting magma--------')
    fig, (ax,ax2,ax3,ax4) = plt.subplots(4,1,figsize=(15,15))
    for kk,model in enumerate(model_list):
        name = 'magma_for_'+model+'.txt'
        filepath = '/home/jiching/geoflac/data/'
        temp1 = np.loadtxt(filepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.plot(yymelt,color=newcolors[kk],label=model)
        ax2.plot(yychamber,color=newcolors[kk],label=model)
#        ax3.bar(fl.time,melt,width=0.1,color=newcolors[kk],label='fmelt')
#        ax4.bar(fl.time,chamber,width=0.1,color=newcolors[kk],label='magma')
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('melt * area',fontsize=20)
    ax2.set_ylabel('chamber *area',fontsize=20)
    ax3.set_ylabel('max melt',fontsize=20)
    ax4.set_ylabel('max magma fraction',fontsize=20)
#    ax.set_xlim(0,time[-1]);ax.grid()
#    ax2.set_xlim(0,time[-1]);ax2.grid()
#    ax3.set_xlim(0,time[-1]);ax3.grid()
#    ax4.set_xlim(0,time[-1]);ax4.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax2.tick_params(axis='x', labelsize=16 )
    ax3.tick_params(axis='x', labelsize=16 )
    ax4.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax2.tick_params(axis='y', labelsize=16 )
    ax3.tick_params(axis='y', labelsize=16 )
    ax4.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    ax.legend()
    fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_magma.png')
    print('=========== DONE =============')
