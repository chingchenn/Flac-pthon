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

fig1=0  ## 4 figure, melting x location + melting phase bar + flat slab length + flat slab depth
fig2=0
fig3=0  ## 2 figure, melting x location + melting phase bar
fig4=1  ## 2 figure, + flat slab length + flat slab depth
fig5=0  ## 3 figure, combine figure
fig6=0  ## 1 figure, melting x location
fig7=0
fig8=1  ## 1 figure, dip only
fig9=0  ## 1 figure, torque only
plt.rcParams["font.family"] = "Times New Roman"
model_list=['Nazca_a0703','Nazca_a0701','Nazca_a0702','Ref_Nazca','Nazca_a0704','Nazca_a0705',]
model_list=['Ref_Nazca']
label_list=['110 km','120 km','130 km','140 km','150 km','160 km']
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
#savepath = '/Users/ji-chingchen/Desktop/data/'
#savepath='D:/model/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
bwith = 3
end = 50
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
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        ax1.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
        ax1.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments+basalt')
    #================================figure setting================================
    # ax2.set_title(model,fontsize=26)
    for aaa in [ax1,ax2,ax3,ax4]:
        
        aaa.tick_params(axis='x', labelsize=16)
        aaa.tick_params(axis='y', labelsize=16)    
        aaa.set_xlim(0,end)
        aaa.grid()
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)

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

    ax5.set_title(model,fontsize=30)
    
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
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
        # divider = make_axes_locatable(ax2)
        # cax = divider.new_vertical(size = '5%', pad = 0.5)
        # cbar=fig.colorbar(qqq,ax=ax2,cax=cax,orientation='horizontal')
        # cbar.set_label('Melting Percentages',fontsize=20)
        # cbar.ax.tick_params(axis='x', labelsize=16)
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        ax1.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='peridotite')
        ax1.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediment+basalt')
    #================================figure setting================================
    # ax2.set_title(model,fontsize=26)
    for aaa in [ax2,ax1]:
        aaa.tick_params(axis='x', labelsize=16)
        aaa.tick_params(axis='y', labelsize=16)
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.grid()
        aaa.set_xlim(0,end)
        
    # ax2.set_ylim(0,400)
    # ax4.set_ylim(-50,-30)
    # ax2.set_ylabel('Distance (km)',fontsize=20)
    # ax1.set_ylabel('molten rocks (km3/km)',fontsize=20)
    # ax1.axes.xaxis.set_visible(False)

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
    for aaa in [ax3,ax4]:
        aaa.tick_params(axis='x', labelsize=16)
        aaa.tick_params(axis='y', labelsize=16)
        aaa.grid()
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.set_xlim(0,end)

    ax3.set_ylim(50,max(length)+20)
    # ax4.set_ylim(50,20)
    # ax4.set_ylim(-50,-30)
    ax3.set_ylabel('Length (km)',fontsize=20)
    ax4.set_ylabel('Depth (km)',fontsize=20)
    ax4.set_xlabel('Time (Myr)',fontsize=20)


if fig5: 
    fig5, (ax)= plt.subplots(3,1,figsize=(10,18))  
    for kk,model in enumerate(model_list):
        ## PLOT Torque 
        time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
        sb = fd.moving_window_smooth(fsb,8)/1e19
        tt = fd.moving_window_smooth(fsu,8)/1e19
        ax[0].plot(time,sb,c='#c06c84',label='slab pull (N)',lw=4)
        ax[0].plot(time,tt,c="#355c7d",label='suction force (N)',lw=4)
        ax[0].legend(fontsize=16,loc='upper left')
        
        ## PLOT melting bar and melting location
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        width=0.15
        ax[1].bar(time,phase_p4+phase_p9,width=width,color='seagreen',label='peridotite')
        ax[1].bar(time,phase_p10,bottom=phase_p4+phase_p9,width=width,color='tomato',label='sediment')
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        
        
        # name = 'magma_for_'+model+'.txt'
        # temp1 = np.loadtxt(savepath+name)
        # melt,chamber,yymelt,yychamber,rrr = temp1.T
        axx1 = ax[1].twinx()
        qqq=axx1.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.03,vmin=0.0)
        # qqq=axx1.scatter(time[melt>0.005],xmelt[melt>0.005],c='#c06c84',s=20)
        # axx1.plot(time,yychamber,color='#282130',lw=4)
        
        
        
        ## PLOT flat slab length & depth & dip!!!
        
        name=model+'_flatslab_time_len.txt'
        time_flat,length,depth=np.loadtxt(savepath+name).T
        smooth_leng= fd.moving_window_smooth(length, 6)
        smooth_dep= fd.moving_window_smooth(depth, 3)
        axx2 = ax[2].twinx()
        lab1=ax[2].plot(time_flat,smooth_leng,label='length of flat slab',color='#2F4F4F',lw=5)
        # axxx2.plot(time,-smooth_dep,label=model,color=newcolors[kk],lw=5)
        
        name = 'plate_dip_of_'+model
        time,dip = np.loadtxt(savepath+name+'.txt').T
        dip_smth = fd.moving_window_smooth(dip[dip>0], 5)
        lab2=axx2.plot(time[dip>0],dip_smth,label='dip',c='royalblue',lw=4)
        lns = lab2+lab1#+lab3#+lab4#+lab5
        labs = [l.get_label() for l in lns]
        
        
        
        ## Add flat slab period
        ax[0].axvspan(time_flat[0],time_flat[-1],facecolor='#35838D', alpha=0.04)
        ax[1].axvspan(time_flat[0],time_flat[-1],facecolor='#35838D', alpha=0.04)
        
        #================================figure setting================================
        ax[0].set_ylabel('Torque (10$^{19}$ N)',fontsize=20)
        ax[1].set_ylim(0,10)
        ax[1].set_ylabel('molten rocks (km3/km)',fontsize=20)
        axx1.tick_params(axis='y', labelsize=16)
        axx1.set_ylabel('Distance (km)',fontsize=20)
        axx1.set_ylim(0,350)
        ax[2].set_ylabel('flat slab depth (km)',fontsize=20)
        # axx2.set_title('Angle Variation of '+str(model),fontsize=24)
        ax[-1].set_xlabel('Time (Myr)',fontsize=20)
        axx2.set_ylabel('Angel ($^\circ$) from '+str(5)+' to '+str(120)+' depth',fontsize=20)
        axx2.tick_params(axis='y', labelsize=16)
        axx2.set_ylim(0,60)
        ax[2].set_ylim(0,150)
        
        for qq in range(len(ax)):
            ax[qq].legend(fontsize=20,facecolor='white',loc='upper left')
            ax[qq].grid()
            ax[qq].spines['bottom'].set_linewidth(bwith)
            ax[qq].spines['top'].set_linewidth(bwith)
            ax[qq].spines['right'].set_linewidth(bwith)
            ax[qq].spines['left'].set_linewidth(bwith)
            ax[qq].tick_params(axis='x', labelsize=16)
            ax[qq].tick_params(axis='y', labelsize=16)
            ax[qq].set_xlim(0,end)
            ax[qq].set_ylim()
        ax[2].legend(lns, labs, fontsize = 20,loc='lower left')
        # ax3.axvspan(time[0],time[-1],facecolor='#35838D', alpha=0.25)
        # ax4.axvspan(time[0],time[-1],facecolor='#35838D', alpha=0.25)
        
        # time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        # qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
    
    # ax2.set_title(model,fontsize=26)
    # fig5.savefig(figpath+model+'_model_result_time.pdf')


if fig6:
    fig6, (ax2)= plt.subplots(1,1,figsize=(10,4))
    for kk,model in enumerate(model_list):
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax2.axvspan(time[0],30,facecolor='#414F67', alpha=0.15)
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax2.scatter(time[melt>0.005],xmelt[melt>0.005],c=melt[melt>0.005],s=65,cmap='OrRd',vmax=0.05,vmin=0.0)
    #================================figure setting================================
    # ax2.set_title(model,fontsize=26)

    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.set_xlim(0,end)
    ax2.set_ylim(0,400)
    ax2.set_ylabel('Distance (km)',fontsize=20)
  #  ax2.set_xlabel('Time (Myr)',fontsize=20)
    ax2.grid()
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    # fig6.savefig('/home/jiching/geoflac/figure/'+str(model)+'metloc.png')

if fig7:
    fig7, (ax6)= plt.subplots(1,1,figsize=(20,16))
    for kk,model in enumerate(model_list):
        # name=model+'_stack_topography.txt'
        # xx,zz = np.loadtxt(savepath+name).T
        name=model+'_stack_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        
        # ax5.plot(xx,zz,c='gray',lw=4)
        # ax6.plot(xs[:-4],zs[:-4],c='k',lw=4)
        # name=model+'_final_topography.txt'
        # xx,zz = np.loadtxt(savepath+name).T
        name=model+'_final_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        withoplot = (xs>0)#*(zs<-10)
        xx = fd.moving_window_smooth(xs[withoplot], 4)
        zz = fd.moving_window_smooth(zs[withoplot], 4)

        
        ax6.plot(xx,-zz,label=label_list[kk],color=newcolors[kk],lw=4)
    #================================figure setting================================

    ax6.legend(fontsize=20)
    ax6.set_xlim(0,900)
    ax6.set_ylim(150,0)
    ax6.set_aspect('equal')
    ymajor_ticks = np.linspace(0,150,num=4)
    ax6.set_yticks(ymajor_ticks)
   
    
    ax6.tick_params(axis='x', labelsize=30)
    ax6.tick_params(axis='y', labelsize=30)
    # ax5.grid()
    ax6.grid()
    # ax6.set_ylabel('Depth (km)',fontsize=20)
    # ax6.set_xlabel('Distance from Ttench (km)',fontsize=20)

    bwith = 3
    ax6.spines['bottom'].set_linewidth(bwith)
    ax6.spines['top'].set_linewidth(bwith)
    ax6.spines['right'].set_linewidth(bwith)
    ax6.spines['left'].set_linewidth(bwith)
    # ax6.legend(fontsize=16)
    # fig7.savefig('/home/jiching/geoflac/figure/'+model+'_forc_flat_melt.png')

if fig8:
    fig8, (ax2)= plt.subplots(1,1,figsize=(10,4))
    for kk,model in enumerate(model_list):
        name = 'plate_dip_of_'+model
        time,dip = np.loadtxt(savepath+name+'.txt').T
        dip_smth = fd.moving_window_smooth(dip[dip>0], 5)
        lab2=ax2.plot(time[dip>0],dip_smth,label=label_list[kk],c=newcolors[kk],lw=4)
        # lns = lab2+lab1#+lab3#+lab4#+lab5
        # labs = [l.get_label() for l in lns]
        #================================figure setting================================
        # ax2.set_title(model,fontsize=26)
        ax2.legend(fontsize=20)
        ax2.tick_params(axis='x', labelsize=16)
        ax2.tick_params(axis='y', labelsize=16)
        ax2.set_xlim(0,end)
        ax2.set_ylim(0,60)
        ax2.set_ylabel('Dip (km)',fontsize=20)
      #  ax2.set_xlabel('Time (Myr)',fontsize=20)
        ax2.grid()
        ax2.spines['bottom'].set_linewidth(bwith)
        ax2.spines['top'].set_linewidth(bwith)
        ax2.spines['right'].set_linewidth(bwith)
        ax2.spines['left'].set_linewidth(bwith)
        # fig6.savefig('/home/jiching/geoflac/figure/'+str(model)+'metloc.png')
if fig9:
    fig9, (ax)= plt.subplots(1,1,figsize=(10,6))  
    for kk,model in enumerate(model_list):
        ## PLOT Torque 
        time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
        sb = fd.moving_window_smooth(fsb,8)/1e19
        tt = fd.moving_window_smooth(fsu,8)/1e19
        ax.plot(time,sb,c='#c06c84',label='slab pull (N)',lw=4)
        ax.plot(time,tt,c="#355c7d",label='suction force (N)',lw=4)
        
         
        ## Add flat slab period
        name=model+'_flatslab_time_len.txt'
        time_flat,length,depth=np.loadtxt(savepath+name).T
        ax.axvspan(time_flat[0],time_flat[-1],facecolor='#35838D', alpha=0.04)
        # ax[1].axvspan(time_flat[0],time_flat[-1],facecolor='#35838D', alpha=0.04)
        
        #================================figure setting================================
        ax.set_ylabel('Torque (10$^{19}$ N)',fontsize=20)
        ax.legend(fontsize=16,loc='upper left')
        ax.set_xlabel('Time (Myr)',fontsize=20)

        
        ax.legend(fontsize=20,facecolor='white',loc='upper left')
        ax.grid()
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_xlim(0,end)
        ax.set_ylim()
       
    fig9.savefig(figpath+model+'_torque_time.pdf')