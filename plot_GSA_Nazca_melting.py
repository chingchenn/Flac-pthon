#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 15:17:00 2023

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
plt.rcParams["font.family"] = "Helvetica"
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'

fig6 = 1 ## Nazca melting phase and melting location v.s. trench
fig8 = 0 ## Nazca 2D magma fraction 
fig6_save = 0
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
    fig6, (ax2,ax1)= plt.subplots(2,1,figsize=(18,12),gridspec_kw={'height_ratios':[1,0.5]})  
    
    model='Nazca_aa06'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    ####====================== FIG1 melting phase ======================
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
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
    cbar.set_label('degree of melting (%)',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    #================================figure setting================================
    #name='Nazca_a0702'+'_flatslab_time_len.txt'
    #time_flat,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        #aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        #aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_ylim(0,8)    
    ax2.set_ylim(0,600)
   # ax3.set_ylim(0,120)
    #ax3.set_ylabel('flat slab length (km)',fontsize=fontsize)
    ax2.set_ylabel('distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('meling rocks (km$^3$/km)',fontsize=fontsize)
    
    ax1.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    if fig6_save:
        fig6.savefig('/Users/chingchen/Desktop/GSA2023/melting_Mexico_time2.pdf')

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

    fig8, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(20,15),gridspec_kw={'height_ratios':[1,1,0.07]})
    xxx_trench = np.max(x_trench)
    chamber_limit = 1e-3
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma_chamber = fl.read_fmagma(i) 
        melt = fl.read_fmelt(i) * 100
        ax2.scatter(ele_x[magma_chamber>chamber_limit],-ele_z[magma_chamber>chamber_limit],c=time_color[i],cmap =cmap,s=10,vmin=0,vmax=40)
        ax1.scatter(ele_x[melt>0.1],-ele_z[melt>0.1],c=time_color[i],cmap = cmap,s=10,vmin=0,vmax=40)


    for ax in [ax1,ax2]:
        ax.grid()
        ax.set_xlim(300,1000)
        ax.set_ylim(150,0)
        ax.set_ylabel('Depth (km)',fontsize=fontsize)
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
        ax1.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        ax2.plot(xx+xxx_trench,-zz,color=time_color[int(frame)*5-1],lw=3)
        
        
    ax2.set_xlabel('Distance (km)',fontsize=fontsize)
    ax1.set_title('location of partial melting',fontsize=fontsize)
    ax2.set_title('location of magma fraction',fontsize=fontsize)
    norm = mpl.colors.Normalize(vmin=0,vmax=40)
    cc1 = mpl.colorbar.ColorbarBase(ax5,cmap=cmap,norm = norm, orientation='horizontal')
    cc1.set_label(label='Time (Myr)', size=fontsize)
    cc1.ax.tick_params(labelsize=fontsize)

    
    if fig8_save:
        fig8.savefig('/Users/chingchen/Desktop/GSA2023/Nazca_magmma2D.pdf')
