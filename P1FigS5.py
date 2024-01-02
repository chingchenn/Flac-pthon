#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 23:52:33 2023

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
fig5S = 1 ## Nazca melting location v.s. trench
fig6_save = 0
fig5S_save = 0

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
    fig6, (ax1,ax2)= plt.subplots(2,1,figsize=(15,18),gridspec_kw={'height_ratios':[1,1]})  
    
    model='Nazca_aa06'
    #model='Nazca_ab01'
    #model='Cocos_aa02'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    ####====================== FIG1 melting phase ======================
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
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
    #melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma=get_magma(end)
    cax = plt.axes([0.99, 0.86, 0.02, 0.279])
    cbar=fig6.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('melting percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    #================================figure setting================================
    name='Nazca_a0702'+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        #aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax1.set_ylim(0,7)    
    ax2.set_ylim(0,700)
   # ax3.set_ylim(0,120)
    #ax3.set_ylabel('flat slab length (km)',fontsize=fontsize)
    ax2.set_ylabel('distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
    
    ax2.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    if fig6_save:
        fig6.savefig('/Users/chingchen/Desktop/FLAC_Works/2023年會/melting_Mexico_time2.pdf')
if fig5S:
    fig5S, (ax2)= plt.subplots(1,1,figsize=(15,9))  
    
    model='Nazca_a0702'
    #model='Nazca_ab04'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    ####====================== melting location ======================
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-8]
        hhh = (melting>1e-8)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100,s=70,cmap='OrRd',vmax=1,vmin=-0)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100,s=70,cmap='OrRd',vmax=1,vmin=-0)
        hod = ele_x[melting>6e-3]
        hhh = (melting>6e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100,s=70,cmap='OrRd',vmax=1,vmin=-0)

    cax = plt.axes([0.99, 0.115, 0.02, 0.856])
    cbar=fig5S.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('degree of melting (%)',fontsize=fontsize)
    cbar.ax.tick_params(axis='y',labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    # #================================figure setting================================
    name=model+'_flatslab_time_len.txt'
    time_flat,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        aaa.set_xlim(0,40)
        aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
    ax2.set_ylim(0,600)
    ax2.set_ylabel('distance from trench (km)',fontsize=fontsize) 
    ax2.set_xlabel('time (Myr)',fontsize=fontsize)
    fig5S.tight_layout()
    if fig5S_save:
        fig5S.savefig(figpath+'FigS5v4.pdf')