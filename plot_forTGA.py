#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 21:07:45 2022

@author: ji-chingchen
"""
import os,sys
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as fd
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig1=0
fig2=0
fig3=0    # Three figure which plot angle, topography and geometry
fig4=0    # Only melting phase vs time
fig5=0    # metloc vs time & melting phase
fig6=0    # Only metloc vs time
fig7=1
plt.rcParams["font.family"] = "Times New Roman"
model_list=['h0924','ch0913']
model_list=['Ref03','h0401','h0402','h0403']
model_list=['h1502','h1504']
model_list=['ch1514','ch1510','ch1509','ch1505','ch1506']
model_list=['ch0913','ch0918','ch0919','ch0920']
bwith = 3
fontsize=25
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
# newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090',
#              '#AE6378','#282130','#7E9680','#24788F',
#              '#849DAB','#EA5E51','#35838D','#4198B9',
#              '#414F67','#97795D','#6B0D47','#A80359',
#              '#52254F']
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/ji-chingchen/Desktop/data/'
savepath = '/Volumes/SSD500/data/'
#savepath='D:/model/data/'
figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'
label_list = ['3$^\circ$','5$^\circ$','7$^\circ$','9$^\circ$']
label_list = ['Flat slab model','Normal slab model']

if fig1:
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


if fig2:
    fig2, (ax)= plt.subplots(1,1,figsize=(14,4))
    for kk,model in enumerate(model_list):
        name='melting_'+model
        time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
        name = 'magma_for_'+model+'.txt'
        temp1 = np.loadtxt(savepath+name)
        melt,chamber,yymelt,yychamber,rrr = temp1.T
        ax.plot(time,yychamber,color=newcolors[kk],label=label_list[kk],lw=4)
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
    fig2.savefig(figpath+str(model)+'metvolume.png')
    # fig4.savefig(figpath+str(model)+'metphase.pdf') 
    #fig2.savefig('D:\\OneDrive - 國立台灣大學/master03/Seminar/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.pdf')
    # fig2.savefig(figpath+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.png')
    fig2.savefig(figpath+'melting_volume'+model_list[0]+'_'+model_list[-1]+'.pdf')
    print('=========== DONE =============') 


if fig3:
    print('--- start plot dip with time ---')
    fig3, (ax)= plt.subplots(3,1,figsize=(8,18),gridspec_kw={'height_ratios':[1,1,0.5]})
    for kk,model in enumerate(model_list):
        ww = kk+3
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax[0].axvspan(time[0],50,ymin=ww/9,ymax=(ww+0.7)/9,facecolor=newcolors[kk], alpha=0.15)
        
        
        name = 'plate_dip_of_'+model
        df = pd.read_csv(savepath+name+'.csv')
        time = df.time[df.angle>0]
        angle = fd.moving_window_smooth(df.angle[df.angle>0],5)
        ax[0].plot(time,angle,lw=5,label=label_list[kk],color=newcolors[kk])
        
        name=model+'_stack_topography.txt'
        xmean,ztop=np.loadtxt(savepath+name).T
        ax[1].plot(xmean,ztop,lw=5,label=label_list[kk],color=newcolors[kk])
        
        xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
        xx= fd.moving_window_smooth(xmean[xmean>0], 6)
        ztop = fd.moving_window_smooth(ztop[xmean>0], 6)
        xmean=xx
        ax[2].plot(xmean,-ztop,c=newcolors[kk],label=model,lw=5)
    ax[0].set_xlim(0,50)
    ax[0].set_title('Angle Variation',fontsize=fontsize)
    ax[0].set_xlabel('Time (Myr)',fontsize=fontsize)
    ax[0].set_ylabel('Angel ($^\circ$) ',fontsize=fontsize)
    ax[2].set_ylim(150,0)
    ax[2].set_xlim(0,400)
    ax[1].legend(fontsize=20,loc='lower right')
    ax[1].set_xlim(-200,600)
    ax[2].set_aspect('equal')
    # ax[1].set_ylim(-2,5)
    for qq in range(len(ax)):
        ax[qq].grid()
        ax[qq].tick_params(axis='x', labelsize=20)
        ax[qq].tick_params(axis='y', labelsize=20)
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
    # fig3.savefig(figpath+'compare_dip_'+model_list[0]+'_'+model_list[-1]+'.png')
    print('=========== DONE =============')
    # fig3.savefig(figpath+str(model)+'sediment_friction_all.pdf') 


if fig4:
    model = model_list[0]
    fig4, (ax)= plt.subplots(1,1,figsize=(10,4))
    name='melting_'+model
    width=0.20
    time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    ax.bar(time,phase_p4+phase_p9,width=width,color='seagreen',label='peridotite')
    ax.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=width,color='tomato',label='sediment+basalt')

    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    ax.axvspan(time[0],30,ymin=2.5/8,ymax=4/8,facecolor=newcolors[11], alpha=0.3)

    #================================figure setting================================
    ymajor_ticks = np.linspace(8,0,num=5)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(0,8)
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
    # fig4.savefig(figpath+str(model)+'metphase.png')
    # fig4.savefig(figpath+str(model)+'metphase.pdf') 



if fig5:
    fig5, (ax)= plt.subplots(2,1,figsize=(10,12))
    for kk,model in enumerate(model_list):
        if kk ==1:
            qq = 1
        else:
            qq = 6
            name='melting_'+model
            time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
            ax[1].bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='peridotite')
            ax[1].bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediment+basalt')
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax[0].axvspan(time[0],30,ymin=qq/8,ymax=(qq+1)/8,facecolor=newcolors[kk], alpha=0.15)
        
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax[0].scatter(time[melt>0.005],xmelt[melt>0.005],c=newcolors[kk],s=70)
        
        
    #================================figure setting================================
    ymajor_ticks = np.linspace(400,0,num=5)
    ax[0].set_yticks(ymajor_ticks)
    ax[0].set_ylim(0,400)
    ax[1].set_ylim(0,8)
    ax[1].set_xlabel('Time (Myr)',fontsize=30)
    ax[0].set_ylabel('Distance (km)',fontsize=30)
    ax[1].set_ylabel('molten rocks (km3/km)',fontsize=30)
    for qq in range(len(ax)):
        ax[qq].tick_params(axis='x', labelsize=26)
        ax[qq].tick_params(axis='y', labelsize=26)
        ax[qq].set_xlim(0,30)
        ax[qq].grid()

        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
    # ax[0].axes.xaxis.set_visible(False)
    # fig5.savefig(figpath+str(model)+'metbot.png')
    # fig5.savefig(figpath+str(model)+'metbot.pdf')
if fig6:
    fig6, (ax)= plt.subplots(1,1,figsize=(10,6))
    for kk,model in enumerate(model_list):
        if kk ==1:
            qq = 1
        else:
            qq = 6
        name=model+'_flatslab_time_len.txt'
        time,length,depth=np.loadtxt(savepath+name).T
        ax.axvspan(time[0],30,ymin=qq/8,ymax=(qq+1)/8,facecolor=newcolors[kk], alpha=0.15)
        
        time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
        qqq=ax.scatter(time[melt>0.005],xmelt[melt>0.005],c=newcolors[kk],s=70)
    #================================figure setting================================
    ymajor_ticks = np.linspace(400,0,num=5)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(0,400)
    ax.set_xlabel('Time (Myr)',fontsize=30)
    ax.set_ylabel('Distance (km)',fontsize=30)

    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.set_xlim(0,30)
    ax.grid()

    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    # ax[0].axes.xaxis.set_visible(False)
    # fig6.savefig(figpath+str(model)+'metloc.png')
    # fig6.savefig(figpath+str(model)+'metloc.pdf')
if fig7:
    fig7, (ax6)= plt.subplots(1,1,figsize=(20,16))
    for kk,model in enumerate(model_list):
        # name=model+'_stack_topography.txt'
        # xx,zz = np.loadtxt(savepath+name).T
        name=model+'_stack_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        name=model+'_final_slab.txt'
        xs,zs = np.loadtxt(savepath+name).T
        withoplot = (xs>0)#*(zs<-10)
        xx = fd.moving_window_smooth(xs[withoplot], 6)
        zz = fd.moving_window_smooth(zs[withoplot], 6)
        if model=='ch0919' or model=='ch0918':
            print(model)
            mx = xx.tolist()
            mz = zz.tolist()
            mx.append(118)
            mz.append(-170)
            xx = np.array(mx)
            zz = np.array(mz)
        
        ax6.plot(xx,-zz,label=model,color=newcolors[kk],lw=5)
    #================================figure setting================================

    # ax5.set_title(model,fontsize=26)
    ax6.set_xlim(0,600)
    ax6.set_ylim(150,0)
    ax6.set_aspect('equal')
    ymajor_ticks = np.linspace(0,150,num=4)
    ax6.set_yticks(ymajor_ticks)
    # ax5.set_aspect('equal')
    
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
    # fig7.savefig(figpath+model+'_geometry.png')
    # fig7.savefig(figpath+model+'_geometry.pdf')

