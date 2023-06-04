#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 21:58:37 2022

@author: chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
#import gravity as fg
import matplotlib
import matplotlib as mpl
#matplotlib.use('Agg')
from matplotlib import cm
from netCDF4 import Dataset
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
#import flac_interpolate as fi


#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Cocos/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
frame_list=[50,70,90,120]
frame = 55
model = 'Nazca_a0702'
model = 'Ref_Cocos'
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
end = fl.nrec

topo_chile          = 0
top_phase           = 0
arc_mexico          = 0
arc_chile           = 0
melting_mexico      = 0
sxz_chile           = 0
sxx_chile           = 1
sxz_mexico          = 0


if topo_chile:   
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    cmap = plt.cm.get_cmap('gist_earth')
    # zmax, zmin =1.5, -10
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(6,6))
    
    for i in range(end):
        x, z = fl.read_mesh(i+1)
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(xt.shape)
        t[:] = i*0.2
        ax.scatter(xt,t,c=zt*4,cmap=cmap)#,vmin=zmin, vmax=zmax)
        trench_t[i] = t[0]
        trench_x[i] = xt[np.argmin(zt)]
        arc_x[i] = xt[np.argmax(zt)]
        # print(np.min(zt))
    ind_within=(arc_x<900)*(trench_x>100)
    # ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(0,1200)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    # distance=arc_x-trench_x
    #ax2.plot(distance[ind_within],trench_t[ind_within],c='b',lw=4)    
    #ax2.set_ylim(0,t[0])
    #ax2.set_xlabel('Distance between arc and trench (km)')
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap)#,vmin=zmin, vmax=zmax)
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    ax_cbin.set_title('Bathymetry (km)')
# plt.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo.jpg')

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
if top_phase:
    cmap = plt.cm.get_cmap('gist_earth')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(10,12))
    
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        ax.scatter(ele_x[:,0],t,c=phase[:,0],cmap=phase15,vmax=20,vmin=1)
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
    ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(450,850)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)")
    ax.set_xlabel("Distance (km)")
    # distance=arc_x-trench_x
    #ax2.plot(distance[ind_within],trench_t[ind_within],c='b',lw=4)    
    #ax2.set_ylim(0,t[0])
    #ax2.set_xlabel('Distance between arc and trench (km)')
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo.jpg')
    
if arc_mexico:
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    # cmap = plt.cm.get_cmap('gist_earth')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(10,6))
    
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        if True in (phase[:,0]==14):
            for kk in range(len(phase[:,0])):
                if phase[kk,0]==14:
                    ax.scatter(ele_x[kk,0]-trench_x[i-1],i*0.2,c='#c98f49',s=70)
                # elif phase[kk,0]==10 and ele_x[kk,0]>arc_x[i-1] and ele_x[kk,0]>670:
                    # ax.scatter(ele_x[kk,0],i*0.2,c='#c98f49',s=10)
        
    # ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    # ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(0,300)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')
    
if arc_chile:
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    # cmap = plt.cm.get_cmap('gist_earth')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(6,6))
    
    # for i in range(126,127+1):
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        # melting = fl.read_fmelt(i)
        melting = fl.read_fmagma(i)*100
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        # if True in (phase[:,0]==14):
        #     for kk in range(len(phase[:,0])):
        #         if phase[kk,0]==14:
        #             ax.scatter(ele_x[kk,0]-trench_x[i-1],i*0.2,c='#c98f49',s=70)
        #         # elif phase[kk,0]==10 and ele_x[kk,0]>arc_x[i-1] and ele_x[kk,0]>670:
        #             # ax.scatter(ele_x[kk,0],i*0.2,c='#c98f49',s=10)


        mmm = 1.5e-3
        if True in (melting>mmm):
            # print(i)
            for kk in range(len(melting[0,:])):
            # for kk in range(51,54):
                if True in (melting[:,kk]>mmm):
                    yy = len(ele_x[:,0][melting[:,kk]>mmm])
                    tt = np.ones(yy)*i*0.2
                    # print(i,kk)
                    ax.scatter(ele_x[:,0][melting[:,kk]>mmm]-trench_x[i-1],tt,c='orange',s=10)#,cmap='OrRd',vmax=0.5,vmin=-0)
                    # print(np.max(melting[:,kk][melting[:,kk]>mmm]*100))
    # ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    # ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(200,600)
    ax.set_ylim(10,30)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')

if melting_mexico:
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    # cmap = plt.cm.get_cmap('gist_earth')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(10,6))
    
    # for i in range(126,127+1):
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        melting = fl.read_fmelt(i)
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        # if True in (phase[:,0]==14):
        #     for kk in range(len(phase[:,0])):
        #         if phase[kk,0]==14:
        #             ax.scatter(ele_x[kk,0]-trench_x[i-1],i*0.2,c='#c98f49',s=70)
        #         # elif phase[kk,0]==10 and ele_x[kk,0]>arc_x[i-1] and ele_x[kk,0]>670:
        #             # ax.scatter(ele_x[kk,0],i*0.2,c='#c98f49',s=10)

        if True in (melting!=0):
            # print(i)
            for kk in range(len(melting[0,:])):
            # for kk in range(51,54):
                if True in (melting[:,kk]!=0):
                    yy = len(ele_x[:,0][melting[:,kk]!=0])
                    tt = np.ones(yy)*i*0.2
                    # print(i,kk)
                    ax.scatter(ele_x[:,0][melting[:,kk]!=0]-trench_x[i-1],tt,c=melting[:,kk][melting[:,kk]!=0]*100,s=10,cmap='OrRd',vmax=3,vmin=-0)
        
    # ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    # ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(0,300)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')
if sxz_chile:
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    cn= plt.cm.get_cmap('seismic')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(6,6))
    
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        # melting = fl.read_fmelt(i)
        sxz = fl.read_sxz(i)*100
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        ax.scatter(ele_x[:,0],t,c=sxz[:,26],cmap=cn,vmax=200,vmin=-200)

    # ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    # ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(0,1200)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')
    
if sxx_chile:
    model = 'Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    cn= plt.cm.get_cmap('seismic')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(10,8))
    
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        sxx = fl.read_sxx(i)*100
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        sssxxx = 0
        for mm in range(0,10):
            sssxxx += sxx[:,mm]
        ax.scatter(ele_x[:,0],t,c=sssxxx/10,cmap=cn,vmax=300,vmin=-300)
        # print(np.max(sssxxx))

    ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    ax.plot(trench_x[ind_within],trench_t[ind_within],c='#7FFF00',lw='4')
    ax.set_xlim(100,1000)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')
    
if sxz_mexico:
    model = 'Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    cn= plt.cm.get_cmap('seismic')
    trench_x = np.zeros(end)
    trench_t = np.zeros(end)
    arc_x = np.zeros(end)
    fig, (ax) = plt.subplots(1,1,figsize=(6,6))
    
    for i in range(1,end+1):
        phase = fl.read_phase(i)
        x, z = fl.read_mesh(i)
        ele_x,ele_z = flac.elem_coord(x, z)
        # melting = fl.read_fmelt(i)
        sxz = fl.read_sxz(i)*100
        xmax = np.amax(x)
        xmin = np.amin(x)
    
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(len(ele_x[:,0]))
        t[:] = i*0.2
        trench_t[i-1] = t[0]
        trench_x[i-1] = xt[np.argmin(zt)]
        arc_x[i-1] = xt[np.argmax(zt)]
        ax.scatter(ele_x[:,0],t,c=sxz[:,16],cmap=cn,vmax=200,vmin=-200)

    # ind_within=(arc_x<900)*(trench_x>100)
    #ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
    # ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
    ax.set_xlim(400,900)
    ax.set_ylim(0,50)
    # ax.set_title(str(model)+" Bathymetry Evolution")
    ax.set_ylabel("Time (Ma)",fontsize = 22)
    ax.set_xlabel("Distance from Trench(km)",fontsize = 22)
    bwith = 3
    # ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    # ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
    # cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    # ax_cbin.set_title('Bathymetry (km)')
    # plt.savefig(figpath+'Mexico_arc_location.pdf')