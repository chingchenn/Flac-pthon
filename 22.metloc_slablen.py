#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 4 17:26:14 2022

@author: ji-chingchen
"""


import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
fig1=1
fig2=0

# model = sys.argv[1]
model = 'ch0810'
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#6B0D47','#A80359','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/ji-chingchen/Desktop/data/'

if fig1: 
    fig, (ax2,ax3,ax4)= plt.subplots(3,1,figsize=(15,12))   
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    ax2.axvspan(time[0],time[-1],ymin=6/8,ymax=7/8,facecolor='#4682B4', alpha=0.8)
    ax3.plot(time,length,c="#000080",lw=3)
    ax4.plot(time,depth,c="#2F4F4F",lw=3)
    time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    qqq=ax2.scatter(time[melt>0],xmelt[melt>0],c=melt[melt>0],cmap='gist_heat_r',vmin=-0,vmax=0.05,s=45)
    # cbar=fig.colorbar(qqq,ax=ax2)
    # cbar.ax.tick_params(labelsize=16)
    #================================figure setting================================
    ax2.set_title(model,fontsize=26)
    
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    
    
    ax2.set_xlim(0,30)
    ax2.set_ylim(0,400)
    ax3.set_xlim(0,30)
    ax3.set_ylim(50,200)
    ax4.set_xlim(0,30)
    ax4.set_ylim(-100,-80)
    
    ax2.set_ylabel('Distance (km)',fontsize=20)
    ax3.set_ylabel('length (km)',fontsize=20)
    ax4.set_ylabel('Depth (km)',fontsize=20)
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    
    ax2.grid()
    ax3.grid()
    ax4.grid()
    
    bwith = 3
    
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
    # fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/master03/Seminar/my present/ch0810.png')
#====================================figure2===================================
if fig2:
    name=model+'_stack_topography.txt'
    xx,zz = np.loadtxt(savepath+name).T
    name=model+'_stack_slab.txt'
    xs,zs = np.loadtxt(savepath+name).T
    fig1, (ax5,ax6)= plt.subplots(2,1,figsize=(12,10))
    # ax5.plot(xx,zz,c='gray',lw=4)
    # ax6.plot(xs[:-4],zs[:-4],c='k',lw=4)
    name=model+'_final_topography.txt'
    xx,zz = np.loadtxt(savepath+name).T
    name=model+'_final_slab.txt'
    xs,zs = np.loadtxt(savepath+name).T
    ax5.plot(xx,zz,c='k',lw=4)
    ax6.plot(xs[zs<0],zs[zs<0],c='k',lw=4)
    ax5.set_title(model,fontsize=26)
    
    ax5.set_xlim(-200,500)
    ax6.set_xlim(-200,500)
    ax6.set_ylim(-170,0)
    ax6.set_aspect('equal')
    # ax5.set_aspect('equal')
    
    ax5.tick_params(axis='x', labelsize=16)
    ax5.tick_params(axis='y', labelsize=16)
    ax6.tick_params(axis='x', labelsize=16)
    ax6.tick_params(axis='y', labelsize=16)
    ax5.grid()
    ax6.grid()
    ax5.set_ylabel('Height (km)',fontsize=20)
    ax6.set_ylabel('Depth (km)',fontsize=20)
    ax6.set_xlabel('Distance from Ttench (km)',fontsize=20)
    # ax2.axes.xaxis.set_visible(False)
    # ax3.axes.xaxis.set_visible(False)
    # ax4.axes.xaxis.set_visible(False)
    
    ax5.spines['bottom'].set_linewidth(bwith)
    ax5.spines['top'].set_linewidth(bwith)
    ax5.spines['right'].set_linewidth(bwith)
    ax5.spines['left'].set_linewidth(bwith)
    ax6.spines['bottom'].set_linewidth(bwith)
    ax6.spines['top'].set_linewidth(bwith)
    ax6.spines['right'].set_linewidth(bwith)
    ax6.spines['left'].set_linewidth(bwith)
    # 
