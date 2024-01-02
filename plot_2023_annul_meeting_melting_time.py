#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  7 23:48:20 2023

@author: chingchen
"""

import flac
import os,sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import function_for_flac as fd
import matplotlib
from matplotlib import cm

#=========================setting=============================
# model = str(sys.argv[1])
# model = 'Nazca_0502'
# path = '/home/jiching/geoflac/'+model+'/'
# path = '/Users/chingchen/Desktop/model/'+model+'/'
path = '/Users/chingchen/Desktop/model/'
plt.rcParams["font.family"] = "Times New Roman"
model='Ref_Cocos'
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'

fig6 = 1 ## Cocos melting phase and melting location v.s. trench
fig7 = 1 ## Cocos melting pahase enlarge
fig8 = 1 ## Nazca 2D magma
fig6_save = 0
fig7_save = 0
fig8_save = 0
#========================= Function =========================
def get_magma(end):
#=========================Time Series=========================
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    rrr=np.zeros(end)
    depth_magma = np.zeros(end)
#=========================main code===========================
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
bwith=3
fontsize = 30
if fig6:
    fig6, (ax2,ax1)= plt.subplots(2,1,figsize=(12,11),gridspec_kw={'height_ratios':[1,1]})  
    
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
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
    melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma=get_magma(end)
    time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    # qqq=ax2.scatter(time[melt>1e-3],xmelt[melt>1e-3],c=melt[melt>1e-3]*100,s=50,cmap='OrRd',vmax=3,vmin=-0)
    # divider = make_axes_locatable(ax2)
    cax = plt.axes([0.99, 0.57, 0.02, 0.37])
    cbar=fig6.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('Melting percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    
    ppptime = time
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    #================================figure setting================================
    
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
        
    ax2.set_ylim(0,300)
    
    ax2.set_ylabel('Distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    
    ax1.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    if fig6_save:
        fig6.savefig('/Users/chingchen/Desktop/FLAC_Works/2023年會/melting_Mexico_time2.pdf')

if fig7:
    fig7, (ax1)= plt.subplots(1,1,figsize=(15,5))  
    # fig7, (ax4,ax5)= plt.subplots(2,1,figsize=(12,8))  
    
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
    
    ppptime = time
    
    
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    #================================figure setting================================
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax1]:
        aaa.tick_params(axis='x', labelsize=fontsize)
        aaa.tick_params(axis='y', labelsize=fontsize)
        aaa.set_xlim(13,40)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    # ax1.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        
    
    ax1.set_ylim(0,0.75)
    
    
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    ax1.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig7.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    if fig7_save:
        fig7.savefig('/Users/chingchen/Desktop/FLAC_Works/2023年會/melting_Mexico_time_enlarge.pdf')



if fig8: 
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

    fig8, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(15,12))
    xxx_trench = np.max(x_trench)
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) * 100
        melt = fl.read_fmelt(i) * 100
        mm=ax2.scatter(ele_x[magma_chamber>1.5e-3],-ele_z[magma_chamber>1.5e-3],color=time_color[i],zorder=1,s=10)
        # time = fl.time[i]
        qqq=ax1.scatter(ele_x[melt>0.1],-ele_z[melt>0.1],color=time_color[i],s = 10)
    def x2dis(x):
        return x -xxx_trench
    def dis2x(x):
        return x +xxx_trench
    ax3 = ax1.secondary_xaxis('top', functions=(x2dis,dis2x))
    ax4 = ax2.secondary_xaxis('top', functions=(x2dis,dis2x))
    #ax3.set_xlabel('Distance from trench (km)',fontsize=30)
    
    for ax in [ax1,ax2]:
        ax.grid()
        ax.set_xlim(300,1000)
        ax.set_ylim(150,0)
        ax.set_ylabel('Depth (km)',fontsize=26)
        ax.tick_params(labelsize=26 )
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(bwith)
        ax.set_aspect('equal')
    for ax in [ax3,ax4]:
        ax.tick_params( labelsize=26 )
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(bwith)
    
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
    
    ax1.set_xlabel('Distance (km)',fontsize=28)
    ax1.set_title('Location of partial melting',fontsize=28)
    ax2.set_title('Location of magma chamber',fontsize=28)
    cax = plt.axes([0.13, 0.29, 0.8, 0.03])
    cc1=fig8.colorbar(cbtime, ax=ax,cax=cax,orientation='horizontal')
    cc1.set_label(label='Time (Myr)', size=23)
    cc1.ax.tick_params(labelsize=26)
    cc1.ax.yaxis.set_label_position('right')
    fig8.tight_layout()
    if fig8_save:
        fig8.savefig('/Users/chingchen/Desktop/FLAC_Works/2023年會/Nazca_magmma2D.pdf')

