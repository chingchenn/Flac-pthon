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
# matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt

#---------------------------------- DO WHAT -----------------------------------
trench_plot             = 0
dip_plot                = 0
plate_geometry          = 0
force_plot_LR           = 0
force_plot_RF           = 0
vel_plot                = 0
stack_topo_plot         = 0
flat_slab_plot          = 0
magma_plot              = 0
melting_phase           = 0
forces                  = 0
forces_suction          = 0
magma_area              = 0
#------------------------------- COMPARE WHAT ---------------------------------
thermal_Nazca = 1
thermal_Cocos = 0
serpentinite_thickness_Nazca = 0
melting_Nazca = 0

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/home/jiching/geoflac/data/'
#savepath = '/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:\\OneDrive - 國立台灣大學/resarch/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
#---------------------------------- SETTING -----------------------------------
if thermal_Nazca:
    model_list=['Nazca_a0701','Nazca_a0702','Ref_Nazca',]#'Nazca_a0704','Nazca_a0705']
    names=['120 km','130 km','140 km','150 km','160 km']
    xmin,xmax = 0,800
    figname='thermal_Nazca'
if thermal_Cocos:
    model_list=['Ref_Cocos','Cocos_a0807']#,'Cocos_a0909','Cocos_a0804']
    names=['20$^\circ$C','15$^\circ$C' ,'17','16']
    xmin,xmax = 0,350
    figname='thermal_Cocos'
if serpentinite_thickness_Nazca:
    model_list=['Nazca_a0910','Nazca_a0702','Nazca_a0706','Nazca_a0707','Nazca_a0708']#,'Nazca_a0709']
    # model_list=['Nazca_a0702','Nazca_a0710','Nazca_a0711','Nazca_a0712','Nazca_a0713']
    # model_list=['Nazca_a0702','Nazca_a0706']
    names=['1','2','3','4','6','10']
    figname = 'serpentinite_thickness_Nazca'
    xmin,xmax = 0,700
if melting_Nazca:
    model_list=['Nazca_a0634','Nazca_a0637','Nazca_a0640','Nazca_a0643',
                'Nazca_a0636','Nazca_a0635',
                'Nazca_a0639','Nazca_a0638',
                'Nazca_a0642','Nazca_a0641',
                'Nazca_a0644','Nazca_a0645',]
                #'Nazca_a0629','Nazca_a0633','Nazca_a0632','Nazca_a0702']
    xmin,xmax = 0,700
    figname = 'melting_Nazca'

# model_list=['Nazca_a0701','Nazca_a0702','Ref_Nazca','Nazca_a0704','Nazca_a0705']
# model_list=['b0607k','b0608k','b0609k','b0611k','b0614k','b0615k']
#model_list=['Cocos_a0646','Cocos_a0807']


newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
plt.rcParams["font.family"] = "Times New Roman"
bwith = 3
##------------------------------------ plot -----------------------------------
if trench_plot:
    print('--- start plotting the trench and topography with time ---')
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    for kk,model in enumerate(model_list):
        name='trench_for_'+model
        df = pd.read_csv(savepath+name+'.csv')
        dis,time,topo=fd.get_topo()
        ax.plot(df.trench_x[df.trench_x>0],df.time[df.trench_x>0],lw=2,label=model,color=newcolors[kk])
    ax.set_xlim(0,dis[-1][-1])
    ax.set_ylim(0,df.time[-1])
    ax.set_title(str(model)+" Bathymetry Evolution",fontsize=24)
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    # fig.savefig(figpath+'compare_trench_location_'+model_list[0]+'_'+model_list[-1]+'.png')
    print('=========== DONE =============')
if dip_plot:
    print('--- start plot dip with time ---')
    fig, (ax2)= plt.subplots(1,1,figsize=(12,6))
    for kk,model in enumerate(model_list):
        name = 'plate_dip_of_'+model
        time,dip = np.loadtxt(savepath+name+'.txt').T
        ddd = fd.moving_window_smooth(dip[dip>0], 3)
        ax2.plot(time[dip>0],ddd,lw=4,label=model,color=newcolors[kk+9])
    ax2.set_xlim(0,50)
    # ax2.set_ylim(0,40)
    ax2.set_ylim(0,60)
    ax2.set_title('Slab Angle Variation',fontsize=24)
    ax2.set_xlabel('Time (Myr)',fontsize=24)
    ax2.set_ylabel('Angel ($^\circ$) ',fontsize=24)
    ax2.grid()
    ax2.legend(fontsize=16)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.tick_params(axis='x', labelsize=24)
    ax2.tick_params(axis='y', labelsize=24)
    fig.savefig(figpath+'compare_dip_'+figname+'.pdf')
    print(figpath+'compare_dip_'+figname+'.pdf')
    print('=========== DONE =============')
if plate_geometry:
    print('--- start plot geometry ---')
    fig2, (ax2) = plt.subplots(1,1,figsize=(14,8))
    for kk,model in enumerate(model_list):  
        xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab_v1.txt').T
        # xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
        print(savepath+str(model)+'_final_slab.txt')
        xx= fd.moving_window_smooth(xmean[xmean>0], 3)
        ztop = fd.moving_window_smooth(ztop[xmean>0], 3)
        xmean=xx
        ax2.plot(xmean,-ztop,c=newcolors[kk+9],label=model,lw=5)
    #ax2.set_xlim(0,max(xmean)+10)
    ax2.set_title("Slab Top Geometry",fontsize=24)
    ax2.set_ylabel("Depth (km)",fontsize=24)
    ax2.set_xlabel("Distance relative to trench (km)",fontsize=24)
    # ax2.legend(fontsize=16)
    xmajor_ticks = np.linspace(0,200,num=5)
    ax2.set_yticks(xmajor_ticks)
    ax2.set_ylim(150,0)
    ax2.set_xlim(xmin,xmax)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.set_aspect('equal')
    ax2.tick_params(axis='x', labelsize=24)
    ax2.tick_params(axis='y', labelsize=24)
    ax2.grid()
    #fig2.savefig('D:\\OneDrive - 國立台灣大學/master03/Seminar/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.pdf')
    # fig2.savefig(figpath+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.pdf')
    fig2.savefig(figpath+'slab_geometry_'+figname+'.pdf')
    print('=========== DONE =============')
if force_plot_LR:
    print('--- start plot left and right force with time ---')
    fig, (ax)= plt.subplots(2,1,figsize=(12,8))   
    for kk,model in enumerate(model_list):
        filepath = savepath+model+'_forc.txt'
        temp1=np.loadtxt(filepath)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        ffl = fd.moving_window_smooth(forc_l,100)
        ffr = fd.moving_window_smooth(forc_r,100)
        ax[0].plot(time,ffl,lw=3,label=model,color=newcolors[kk])
        ax[1].plot(time,ffr,lw=3,label=model,color=newcolors[kk])
    ax[1].set_title('continental side force',fontsize=16)
    ax[0].set_title('oceanic side force',fontsize=16)
    # ax[0].set_ylim(1.5e12,1.5e13)
    for qq in range(len(ax)):
        ax[qq].tick_params(axis='x', labelsize=16)
        ax[qq].tick_params(axis='y', labelsize=16)
        ax[qq].grid()
        ax[qq].set_xlim(0,time[-1])
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
    # ax[1].legend(fontsize=16)
    # fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_forc.pdf')
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
    # fig2.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_ringforc.png')
    print('=========== DONE =============')
if vel_plot:
    print('--- start plot velocity with time ---')
    fig3, (ax4)= plt.subplots(2,1,figsize=(10,8))
    for kk,model in enumerate(model_list):
        filepath = savepath+model+'_forc.txt'
        temp1=np.loadtxt(filepath)
        nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
        ax4[0].plot(time,vl*31545741325,lw=5,label=model,color=newcolors[kk])
        ax4[1].plot(time,vr*31545741325,lw=5,label=model,color=newcolors[kk])
    ax4[1].set_xlabel('Time (Myr)',fontsize=16)
    ax4[0].set_ylabel('Velocity (mm/yr)',fontsize=16)
    ax4[1].set_ylabel('Velocity (mm/yr)',fontsize=16)
    ax4[0].set_title('oceanic side velocity',fontsize=16)
    for qq in range(len(ax4)):
        ax4[qq].set_xlim(0,time[-1])
        ax4[qq].tick_params(axis='x', labelsize=16)
        ax4[qq].tick_params(axis='y', labelsize=16)
        ax4[qq].grid()
        ax4[qq].spines['bottom'].set_linewidth(bwith)
        ax4[qq].spines['top'].set_linewidth(bwith)
        ax4[qq].spines['right'].set_linewidth(bwith)
        ax4[qq].spines['left'].set_linewidth(bwith)
        ax4[qq].legend(fontsize=26)
    # fig3.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_vel.png')
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
    fig2, (ax) = plt.subplots(2,1,figsize=(10,6))
    for kk,model in enumerate(model_list):
       name=model+'_flatslab_time_len.txt'
       time,length,depth=np.loadtxt(savepath+name).T
       ax[0].plot(time,length,lw=3,label=model,color=newcolors[kk])
       ax[1].plot(time,-depth,lw=3,label=names[kk],color=newcolors[kk])
    ax[0].set_ylabel('Length (km)',fontsize=18)
    ax[1].set_xlabel('Time (Myr)',fontsize=18)
    ax[1].set_ylabel('Depth (km) ',fontsize=18)
    ax[1].legend(fontsize=16)
    ax[1].set_ylim(130,80)
    ax[0].set_title('Flat Slab Length',fontsize=20)
    ax[1].set_title('Flat Slab Depth',fontsize=20)
    for qq in range(len(ax)):
        ax[qq].tick_params(axis='x', labelsize=16)
        ax[qq].tick_params(axis='y', labelsize=16)
        ax[qq].grid()
        ax[qq].set_xlim(0,50)
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
    fig2.tight_layout()
    fig2.savefig(figpath+figname+'_flatslab_length.pdf')
    print('=========== DONE =============')
if magma_plot:
    print('--------plotting magma--------')
    fig, (ax) = plt.subplots(4,1,figsize=(15,15))
    for kk,model in enumerate(model_list):
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax[0].plot(yymelt,color=newcolors[kk],label=model,lw=3)
        ax[1].plot(yychamber,color=newcolors[kk],label=model,lw=3)
        ax[2].bar(time,melt,width=0.1,color=newcolors[kk],label='fmelt')
        ax[3].bar(time,chamber,width=0.1,color=newcolors[kk],label='magma')
    ax[3].set_xlabel('Time (Myr)',fontsize=20)
    ax[0].set_ylabel('melt * area',fontsize=20)
    ax[1].set_ylabel('chamber *area',fontsize=20)
    ax[2].set_ylabel('max melt',fontsize=20)
    ax[3].set_ylabel('max magma fraction',fontsize=20)
    ax[0].set_title('Model : '+model,fontsize=25)
    for qq in range(len(ax)):
        ax[qq].set_xlim(0,time[-1])
        ax[qq].grid()
        ax[qq].tick_params(axis='x', labelsize=16 )
        ax[qq].tick_params(axis='y', labelsize=16 )
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
    ax[0].legend()
    # fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_magma.png')
    print('=========== DONE =============')
if melting_phase:
    fig4, (ax)= plt.subplots(1,1,figsize=(10,4))
    for kk,model in enumerate(model_list):
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.plot(time,yychamber,color=newcolors[kk],label=model)
#    for kk,model in enumerate(model_list):
#        name='melting_'+model
#        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
#        total = np.zeros(len(time))
#        total_volume = phase_p4+phase_p3+phase_p10
#        ax.plot(time,total_volume,lw=3,label=model,color=newcolors[kk])
    #ymajor_ticks = np.linspace(8,0,num=5)
    #ax.set_yticks(ymajor_ticks)
    #ax.set_ylim(0,8)
    ax.set_xlabel('Time (Myr)',fontsize=30)
    # ax.set_ylabel('molten rocks (km3/km)',fontsize=30)

    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.set_xlim(0,30)
    ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.legend(fontsize = 20)
    # fig4.savefig(figpath+str(model)+'metvolume.png')
    # fig4.savefig(figpath+str(model)+'metphase.pdf') 
if forces:
    fig, (ax)= plt.subplots(1,1,figsize=(12,8))
    ## PLOT Torque 
    for kk,model in enumerate(model_list):
        time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
        sb = fd.moving_window_smooth(fsb,8)/1e18
        tt = fd.moving_window_smooth(fsu,8)/1e18
        ax.plot(time,sb,c=newcolors[kk],label='slab pull (N)',lw=3)
        ax.plot(time,tt,c=newcolors[kk],label='suction force (N)',lw=3)

    #================================figure setting================================
    ax.set_ylabel('Torque (10$^{18}$ N)',fontsize=20)
    # ax.legend(fontsize=16,loc='upper left')
    ax.set_xlabel('Time (Myr)',fontsize=20)

    # ax.legend(fontsize=20,facecolor='white',loc='upper left')
    ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlim(0,30)
    ax.set_ylim()
    
if forces_suction:
    fontsize=30
    fig, (ax)= plt.subplots(1,1,figsize=(15,7))
    ## PLOT Torque 
    for kk,model in enumerate(model_list):
        time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
        # sb = fd.moving_window_smooth(fsb,8)/1e18
        tt = fd.moving_window_smooth(fsu,6)/1e19
        # ax.plot(time,sb,c=newcolors[kk+1],label='slab pull (N)',lw=3,ls='--')
        ax.plot(time,tt,c=newcolors[kk+9],label=names[kk],lw=4)

    #================================figure setting================================
    ax.set_ylabel('Torque (10$^{19}$ N)',fontsize=fontsize)
    # ax.legend(fontsize=16,loc='upper left')
    ax.set_xlabel('Time (Myr)',fontsize=fontsize)

    # ax.legend(fontsize=fontsize,facecolor='white',loc='upper left')
    ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    # ax.set_xlim(0,40)
    if serpentinite_thickness_Nazca:
        ax.set_ylim(-0.5,3)
        ax.set_xlim(0,40)
    if thermal_Cocos:
        ax.set_ylim(-0.5,1.5)
        ax.set_xlim(0,30)
    if thermal_Nazca:
        ax.set_ylim(-1,9)
        ax.set_xlim(0,40)
    fig.savefig(figpath+figname+'_suction.pdf')
if magma_area:
    print('--------plotting magma--------')
    fig, (ax) = plt.subplots(1,1,figsize=(12,8))
    for kk,model in enumerate(model_list):
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.bar(time,yychamber,width=0.2,color=newcolors[kk],label=model,lw=3)
    ax.set_ylabel('km$^3$/km',fontsize=30)
    ax.set_xlim(0,30)
    ax.grid()
    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.set_xlabel('Time (Myr)',fontsize=30)
    # ax.legend(fontsize=24)
    # fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_magma.png')
    print('=========== DONE =============')
if melting_Nazca:
    fontsize = 30
    fig, (ax2,ax)= plt.subplots(2,1,figsize=(16,13),gridspec_kw={'height_ratios':[0.5,1]})
    for kk,model in enumerate(model_list):
        
        if kk <=3:
            ccc = newcolors[1]
            qqq = 'not flat slab'
        else:
            ccc = newcolors[4]
            qqq = 'flat slab'
        
        xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
        xx= fd.moving_window_smooth(xmean[xmean>0], 5)
        ztop = fd.moving_window_smooth(ztop[xmean>0], 5)
        xmean=xx
        ax2.plot(xmean,-ztop,c=ccc,label=qqq,lw=5)
        
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.plot(time,yychamber,color=ccc,label=qqq,lw=3)
        
    for nnn in [ax,ax2]:
        nnn.spines['bottom'].set_linewidth(bwith)
        nnn.spines['top'].set_linewidth(bwith)
        nnn.spines['right'].set_linewidth(bwith)
        nnn.spines['left'].set_linewidth(bwith)
        nnn.tick_params(axis='x', labelsize=fontsize)
        nnn.tick_params(axis='y', labelsize=fontsize)
        nnn.grid()
    ax2.set_title('Slab Geometry',fontsize=fontsize+10)
    ax2.set_xlabel("Distance relative to trench (km)",fontsize=fontsize)
    ax2.set_ylabel("Depth (km)",fontsize=fontsize)
    ax.set_xlabel('Time (Myr)',fontsize=fontsize)
    ax.set_ylabel('Volume (km$^3$/km)',fontsize=fontsize+5)
    ax.set_title('Magma chamber volume v.s. Time',fontsize=fontsize+10)
    xmajor_ticks = np.linspace(0,200,num=5)
    ax2.set_yticks(xmajor_ticks)
    ax2.set_ylim(150,0)
    ax2.set_xlim(xmin,xmax)
    ax2.set_aspect('equal')
    ax.set_xlim(0,30)
    ax.set_ylim(0,50)
    
    fig.tight_layout()
    # fig.savefig(figpath+'magma_area_compare.pdf')
    figwww, (ax2,ax,ax4)= plt.subplots(3,1,figsize=(16,20),gridspec_kw={'height_ratios':[0.5,1,1]})
    for kk,model in enumerate(model_list):
        
        if kk <=3:
            ccc = newcolors[3]
            qqq = 'not flat slab'
        elif kk > 7:
            ccc = newcolors[4]
            qqq = 'flat slab'
        else:
            ccc = newcolors[5]
            qqq = 'flat slab'
        
        xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
        xx= fd.moving_window_smooth(xmean[xmean>0], 5)
        ztop = fd.moving_window_smooth(ztop[xmean>0], 5)
        xmean=xx
        ax2.plot(xmean,-ztop,c=ccc,label=qqq,lw=5)
        
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.plot(time,yychamber,color=ccc,label=qqq,lw=3)
        
        # name = 'plate_dip_of_'+model
        # time,dip = np.loadtxt(savepath+name+'.txt').T
        # ddd = fd.moving_window_smooth(dip[dip>0], 3)
        # ax3.plot(time[dip>0],ddd,lw=4,label=qqq,color=ccc)
        
        
        # PLOT Torque 
        
        time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
        # sb = fd.moving_window_smooth(fsb,8)/1e18
        tt = fd.moving_window_smooth(fsu,6)/1e19
        # ax.plot(time,sb,c=newcolors[kk+1],label='slab pull (N)',lw=3,ls='--')
        ax4.plot(time,tt,c=ccc,label=qqq,lw=4)

        #================================figure setting================================
        ax4.set_ylabel('Torque (10$^{19}$ N)',fontsize=fontsize)
        # ax.legend(fontsize=16,loc='upper left')
        ax4.set_xlabel('Time (Myr)',fontsize=fontsize)
        
    for nnn in [ax,ax2,ax4]:
        nnn.spines['bottom'].set_linewidth(bwith)
        nnn.spines['top'].set_linewidth(bwith)
        nnn.spines['right'].set_linewidth(bwith)
        nnn.spines['left'].set_linewidth(bwith)
        nnn.tick_params(axis='x', labelsize=fontsize)
        nnn.tick_params(axis='y', labelsize=fontsize)
        nnn.grid()
        
    # ax3.set_xlim(0,30)
    # ax3.set_ylim(0,60)
    # ax3.set_title('Slab Angle Variation',fontsize=fontsize)
    # ax3.set_xlabel('Time (Myr)',fontsize=fontsize)
    # ax3.set_ylabel('Angel ($^\circ$) ',fontsize=fontsize)
    ax2.set_title('Slab Geometry',fontsize=fontsize+10)
    ax2.set_xlabel("Distance relative to trench (km)",fontsize=fontsize)
    ax2.set_ylabel("Depth (km)",fontsize=fontsize)
    ax.set_xlabel('Time (Myr)',fontsize=fontsize)
    ax.set_ylabel('Volume (km$^3$/km)',fontsize=fontsize+5)
    ax.set_title('Magma chamber volume v.s. Time',fontsize=fontsize+10)
    xmajor_ticks = np.linspace(0,200,num=5)
    ax2.set_yticks(xmajor_ticks)
    ax2.set_ylim(150,0)
    ax2.set_xlim(xmin,xmax)
    ax2.set_aspect('equal')
    ax.set_xlim(0,30)
    ax4.set_xlim(0,30)
    
    figwww.tight_layout()
    # figwww.savefig(figpath+'magma_area_compare.pdf')
