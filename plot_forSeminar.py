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
import function_for_flac as fd
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig1=0
fig2=0
fig3=1
fig4=1
plt.rcParams["font.family"] = "Times New Roman"
model_list=['h0829']

newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/ji-chingchen/Desktop/data/'
#savepath='D:/model/data/'


if fig1: 
    fig, (ax2,ax1,ax3,ax4)= plt.subplots(4,1,figsize=(14,10))  
    for kk,model in enumerate(model_list):
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax1.axvspan(time[0],time[-1],facecolor='#A80359', alpha=0.2)
        ax2.axvspan(time[0],time[-1],facecolor='#35838D', alpha=0.25)
        ax3.axvspan(time[0],time[-1],facecolor='#35838D', alpha=0.25)
        ax4.axvspan(time[0],time[-1],facecolor='#35838D', alpha=0.25)
        smooth_leng= fd.moving_window_smooth(length, 6)
        smooth_dep= fd.moving_window_smooth(depth, 3)
        ax3.plot(time,smooth_leng,label=model,color=newcolors[kk],lw=5)
        ax4.plot(time,smooth_dep,label=model,color=newcolors[kk],lw=5)
        # ax3.plot(time,length,label=model,color=newcolors[kk],lw=3)
        # ax4.plot(time,depth,label=model,color=newcolors[kk],lw=3)
        time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        ax1.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
        ax1.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments+basalt')
    #================================figure setting================================
    # ax2.set_title(model,fontsize=26)
    
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)

    
    ax4.set_xlim(0,30)
    ax1.set_xlim(0,30)
    ax2.set_xlim(0,30)
    ax3.set_xlim(0,30)
    ax2.set_ylim(0,400)
    ax3.set_ylim(50,max(length)+20)
    ax4.set_ylim(np.average(depth)-10,np.average(depth)+10)
    # ax4.set_ylim(-50,-30)
    
    ax2.set_ylabel('Distance (km)',fontsize=20)
    ax3.set_ylabel('Length (km)',fontsize=20)
    ax4.set_ylabel('Depth (km)',fontsize=20)
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    ax1.set_ylabel('molten rocks (km3/km)',fontsize=20)
    # ax1.axes.xaxis.set_visible(False)
    
    ax1.grid()
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
    ax1.spines['bottom'].set_linewidth(bwith)
    ax1.spines['top'].set_linewidth(bwith)
    ax1.spines['right'].set_linewidth(bwith)
    ax1.spines['left'].set_linewidth(bwith)
    ax1.legend(fontsize=25)

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



if fig3: 
    fig, (ax2,ax1)= plt.subplots(2,1,figsize=(10,10),gridspec_kw={'height_ratios':[1,1]})  
    for kk,model in enumerate(model_list):
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax1.axvspan(time[0],time[-1],facecolor='#414F67', alpha=0.2)
        ax2.axvspan(time[0],time[-1],facecolor='#414F67', alpha=0.15)
        time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
        # divider = make_axes_locatable(ax2)
        # cax = divider.new_vertical(size = '5%', pad = 0.5)
        # cbar=fig.colorbar(qqq,ax=ax2,cax=cax,orientation='horizontal')
        # cbar.set_label('Melting Percentages',fontsize=20)
        # cbar.ax.tick_params(axis='x', labelsize=16)
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        ax1.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
        ax1.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments+basalt')
    #================================figure setting================================
    # ax2.set_title(model,fontsize=26)
    
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)

    
    ax1.set_xlim(0,30)
    ax2.set_xlim(0,30)
    ax2.set_ylim(0,400)

    # ax4.set_ylim(-50,-30)
    
    # ax2.set_ylabel('Distance (km)',fontsize=20)
    # ax1.set_ylabel('molten rocks (km3/km)',fontsize=20)
    # ax1.axes.xaxis.set_visible(False)
    
    ax1.grid()
    ax2.grid()

    
    bwith = 3

    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax1.spines['bottom'].set_linewidth(bwith)
    ax1.spines['top'].set_linewidth(bwith)
    ax1.spines['right'].set_linewidth(bwith)
    ax1.spines['left'].set_linewidth(bwith)
    ax1.legend(fontsize=22,facecolor='white')

if fig4: 
    fig3, (ax3,ax4)= plt.subplots(2,1,figsize=(10,10))  
    for kk,model in enumerate(model_list):
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        smooth_leng= fd.moving_window_smooth(length, 6)
        smooth_dep= fd.moving_window_smooth(depth, 3)
        ax3.plot(time,smooth_leng,label=model,color=newcolors[kk],lw=5)
        ax4.plot(time,-smooth_dep,label=model,color=newcolors[kk],lw=5)

    #================================figure setting================================
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)


    
    ax4.set_xlim(0,30)
    ax3.set_xlim(0,30)

    ax3.set_ylim(50,max(length)+20)
    ax4.set_ylim(50,20)
    # ax4.set_ylim(-50,-30)
    

    # ax3.set_ylabel('Length (km)',fontsize=20)
    # ax4.set_ylabel('Depth (km)',fontsize=20)
    # ax4.set_xlabel('Time (Myr)',fontsize=20)


    ax3.grid()
    ax4.grid()
    
    bwith = 3
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)

