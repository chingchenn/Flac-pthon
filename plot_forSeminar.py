#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 13:46:35 2022

@author: ji-chingchen
"""

import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
fig1=0
fig2=1

model_list=['ch0810','ch0811']
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/ji-chingchen/Desktop/data/'


if fig1: 
    fig, (ax2,ax3,ax4)= plt.subplots(3,1,figsize=(15,12))  
    for kk,model in enumerate(model_list):
        name = model+'_forc.txt'
        temp1=np.loadtxt(savepath+name)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        # ax.scatter(time,forc_l,label=model,color='#2F4F4F',s=5)
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        # ax2.axvspan(time[0],time[-1],facecolor='#52254F', alpha=0.1)
        ax3.plot(time,length,label=model,color=newcolors[kk],lw=3)
        ax4.plot(time,depth,label=model,color=newcolors[kk],lw=3)
        time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        ax2.scatter(time[melt>0],xmelt[melt>0],label=model,color=newcolors[kk],s=45)
        
    #================================figure setting================================
    ax2.set_title(model,fontsize=26)
    
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    
    ax4.set_xlim(0,30)
    ax2.set_xlim(0,30)
    ax3.set_xlim(0,30)
    ax2.set_ylim(0,400)
    ax3.set_ylim(50,200)
    ax4.set_ylim(-100,-80)
    
    ax2.set_ylabel('Distance (km)',fontsize=20)
    ax3.set_ylabel('length (km)',fontsize=20)
    ax4.set_ylabel('Depth (km)',fontsize=20)
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    
    ax2.grid()
    ax3.grid()
    ax4.grid()
    
    bwith = 3
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax4.legend(fontsize=16)

#====================================figure2===================================
if fig2:
    fig1, (ax5,ax6)= plt.subplots(2,1,figsize=(15,12))
    for kk,model in enumerate(model_list):
        name=model+'_stack_topography.txt'
        xx,zz = np.loadtxt(savepath+name).T
        name=model+'_stack_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        
        # ax5.plot(xx,zz,c='gray',lw=4)
        # ax6.plot(xs[:-4],zs[:-4],c='k',lw=4)
        name=model+'_final_topography.txt'
        xx,zz = np.loadtxt(savepath+name).T
        name=model+'_final_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        ax5.plot(xx,zz,label=model,color=newcolors[kk],lw=4)
        ax6.plot(xs[zs<0],zs[zs<0],label=model,color=newcolors[kk],lw=4)
    #================================figure setting================================

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
    bwith = 3
    ax5.spines['bottom'].set_linewidth(bwith)
    ax5.spines['top'].set_linewidth(bwith)
    ax5.spines['right'].set_linewidth(bwith)
    ax5.spines['left'].set_linewidth(bwith)
    ax6.spines['bottom'].set_linewidth(bwith)
    ax6.spines['top'].set_linewidth(bwith)
    ax6.spines['right'].set_linewidth(bwith)
    ax6.spines['left'].set_linewidth(bwith)
    ax6.legend(fontsize=16)
    # fig.savefig('/home/jiching/geoflac/figure/'+model+'_forc_flat_melt.png')
