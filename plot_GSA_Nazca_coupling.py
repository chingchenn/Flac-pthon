#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:30:14 2023

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



plt.rcParams["font.family"] = "Arial"
#---------------------------------- DO WHAT -----------------------------------

### pdf or png
png             = 0
pdf             = 0

### plot
plot_topo_and_xvel=0
plot_topo_and_zvel=0
plot_topo = 1
plot_sxx = 0
plot_sII = 0
plot_s1_profile = 0
plot_strength_litho = 0
mp4=0

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
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Nazca/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
# figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/AGU/POSTER/Poster_figure/'



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

xmin,xmax=250,1000
zmin,zmax=-200,10
model = 'Nazca_aa06'
shift = 320
#model = 'Ref_Cocos'
frame = 51
frame = 11
frame = 190
frame = 126
frame = 200
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

    
g=10
def compute_s1(sxx, szz, sxz):
    mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
    theta = 0.5 * np.arctan2(2*sxz,  sxx-szz)

    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    nx, nz = sxx.shape
    tmp = np.zeros((nx, nz, 3), dtype=sxx.dtype)
    tmp[:,:,0] = mag * np.sin(theta)
    tmp[:,:,1] = mag * np.cos(theta)
    return tmp


fig, (ax)= plt.subplots(1,1,figsize=(17,10))
temp = fl.read_temperature(frame)
bwith = 3
#--------------------- plotting -------------------------
x,z = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(x, z)
ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
sxx = fl.read_sxx(frame)*100
cbsxx = plt.cm.get_cmap('seismic_r')
cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
ax.set_title('sxx',fontsize=25)
cax = plt.axes([0.945, 0.365, 0.01, 0.271])
cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
cc1.ax.tick_params(labelsize=20)
ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
cc1.ax.yaxis.set_label_position('left')
ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
#ax.scatter(ele_x[sxx>50],-ele_z[sxx>50],c='green')
# ---------------------- plot setting --------------------------
ax.set_aspect('equal')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(bwith)
ax.tick_params(axis='x', labelsize=23)
ax.tick_params(axis='y', labelsize=23)
ymajor_ticks = np.linspace(300,0,num=7)
ax.set_yticks(ymajor_ticks)
ax.set_ylim(-zmin,-zmax)
ax.set_xlim(xmin,xmax)
ax.set_xlabel('distance (km)',fontsize=26)
ax.set_ylabel('depth (km)',fontsize=26)

ax.set_xlim(xmin,xmax)
# ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
xx=xx[zz<0]
zz=zz[zz<0]
zz = fd.moving_window_smooth(zz,3)
ax.plot(xx+shift,-zz,color='k',lw=5)

ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
#fig.savefig(figpath+model+'frame_'+str(frame)+'_sxx_gourand.pdf')


end = 228
if plot_topo_and_xvel:
    for frame in range(1,end):
    
        fig2, (ax2,ax)= plt.subplots(2,1,figsize=(15,10),gridspec_kw={'height_ratios':[2,3]})
      
        ax3 = ax2.twinx()
        ax3.plot(x[:,0],z[:,0]+2,lw = 3,color = '#4198B9')
        
        xvel,zvel= fl.read_vel(frame)
        cbsxx = plt.cm.get_cmap('seismic')
        smooth_vel = fd.moving_window_smooth(xvel[:,0], 3)
        ax2.plot(x[:,0],smooth_vel,lw =3,color = '#455E45')
        
        
        temp = fl.read_temperature(frame)
        
        #--------------------- plotting -------------------------
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax.set_title('sxx',fontsize=25)
        cax = plt.axes([0.945, 0.365, 0.01, 0.271])
        cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
        cc1.ax.tick_params(labelsize=20)
        ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        # ---------------------- plot setting --------------------------
        ax.set_aspect('equal')
        for aa in [ax3,ax2,ax]:
            aa.set_xlim(xmin,xmax)
            aa.tick_params(labelsize=23)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
            
        ymajor_ticks = np.linspace(300,0,num=7)
        ax.set_yticks(ymajor_ticks)
        ax.set_ylim(-zmin,-zmax)
        ax2.set_ylim(-0.1,0.2)
        ax3.set_ylim(-2,4)
        ax3.tick_params(labelcolor="#4198B9")
        ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
        
        if frame < 10:
            qq = '00'+str(frame)
        elif frame < 100 and frame >=10:
            qq = '0'+str(frame)
        else:
            qq=str(frame)
        fig2.savefig('/Users/chingchen/Desktop/GSA2023/mp4file/frame_'+qq+'_Vx.png')
    if mp4: 
        from PIL import Image
        import glob
         
        # Create the frames
        frames = []
        for i in  range(1,end):
            if i < 10:
                qq = '00'+str(i)
            elif i < 100 and i >=10:
                qq = '0'+str(i)
            else:
                qq=str(i)
            
            img='/Users/chingchen/Desktop/GSA2023/mp4file/frame_'+qq+'_Vx.png'
            new_frame = Image.open(img)
            frames.append(new_frame)
         
        # Save into a GIF file that loops forever
        frames[0].save('/Users/chingchen/Desktop/GSA2023/mp4file/snapshot_Vx.gif', format='GIF', append_images=frames[1:], 
                       save_all=True, duration=80, loop=0)
    #-----------------------------creat mp4-----------------------------------------    
        import moviepy.editor as mp
        clip = mp.VideoFileClip('/Users/chingchen/Desktop/GSA2023/mp4file/snapshot_Vx.gif')
        #clip.write_videofile(figpath+'phase_vis_'+model+".mp4")
        clip.write_videofile('/Users/chingchen/Desktop/GSA2023/mp4file/snapshotVx.mp4')

if plot_topo_and_zvel:
    for frame in range(1,end,2):
    
        fig3, (ax2,ax)= plt.subplots(2,1,figsize=(15,10),gridspec_kw={'height_ratios':[2,3]})
      
        ax3 = ax2.twinx()
        ax3.plot(x[:,0],z[:,0]+2,lw = 3,color = '#4198B9')
        
        xvel,zvel= fl.read_vel(frame)
        cbsxx = plt.cm.get_cmap('seismic')
        smooth_vel = fd.moving_window_smooth(zvel[:,0], 3)
        ax2.plot(x[:,0],smooth_vel,lw =3,color = '#455E45')
        
        
        temp = fl.read_temperature(frame)
        
        #--------------------- plotting -------------------------
        x,z = fl.read_mesh(frame)
        ele_x,ele_z = flac.elem_coord(x, z)
        ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
        sxx = fl.read_sxx(frame)*100
        cbsxx = plt.cm.get_cmap('seismic')
        cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx,cmap=cbsxx,vmin=-300,vmax=300,shading='gouraud')
        ax.set_title('sxx',fontsize=25)
        cax = plt.axes([0.945, 0.365, 0.01, 0.271])
        cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
        cc1.ax.tick_params(labelsize=20)
        ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
        cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
        cc1.ax.yaxis.set_label_position('left')
        ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        # ---------------------- plot setting --------------------------
        ax.set_aspect('equal')
        for aa in [ax3,ax2,ax]:
            aa.set_xlim(xmin,xmax)
            aa.tick_params(labelsize=23)
            for axis in ['top','bottom','left','right']:
                aa.spines[axis].set_linewidth(bwith)
            
        ymajor_ticks = np.linspace(300,0,num=7)
        ax.set_yticks(ymajor_ticks)
        ax.set_ylim(-zmin,-zmax)
        ax2.set_ylim(-0.1,0.2)
        ax3.set_ylim(-2,4)
        ax3.tick_params(labelcolor="#4198B9")
        ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
        
        if frame < 10:
            qq = '00'+str(frame)
        elif frame < 100 and frame >=10:
            qq = '0'+str(frame)
        else:
            qq=str(frame)
        fig3.savefig('/Users/chingchen/Desktop/GSA2023/mp4file/frame_'+qq+'_Vz.png')
    if mp4: 
        from PIL import Image
        import glob
         
        # Create the frames
        frames = []
        for i in  range(1,end):
            if i < 10:
                qq = '00'+str(i)
            elif i < 100 and i >=10:
                qq = '0'+str(i)
            else:
                qq=str(i)
            
            img='/Users/chingchen/Desktop/GSA2023/mp4file/frame_'+qq+'_Vz.png'
            new_frame = Image.open(img)
            frames.append(new_frame)
         
        # Save into a GIF file that loops forever
        frames[0].save('/Users/chingchen/Desktop/GSA2023/mp4file/snapshot_Vz.gif', format='GIF', append_images=frames[1:], 
                       save_all=True, duration=80, loop=0)
    #-----------------------------creat mp4-----------------------------------------    
        import moviepy.editor as mp
        clip = mp.VideoFileClip('/Users/chingchen/Desktop/GSA2023/mp4file/snapshot_Vz.gif')
        #clip.write_videofile(figpath+'phase_vis_'+model+".mp4")
        clip.write_videofile('/Users/chingchen/Desktop/GSA2023/mp4file/snapshot.mp4')



newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 


if plot_topo:
    cmap = plt.cm.get_cmap('gist_earth')
    fig4, (ax)= plt.subplots(1,1,figsize=(8,10))
    tip_x = np.zeros(end)
    tip_z = np.zeros(end)
    tttime = np.zeros(end)
    x_change_ec = np.zeros(end)
    z_change_ec = np.zeros(end)
    x_change_serp = np.zeros(end)
    z_change_serp = np.zeros(end)
    x_change_cho = np.zeros(end)
    z_change_cho = np.zeros(end)
    for frame in range(1,end):
        x, z = fl.read_mesh(frame)
        ele_x, ele_z = flac.elem_coord(x,z)
        magma = fl.read_fmelt(frame)
        phase = fl.read_phase(frame)
        xt = x[:,0]
        t = np.zeros(xt.shape)
        t[:] = frame*0.2
        melt_thod = 1e-3
        sxx = fl.read_sxx(frame)*100
        cb_plot=ax.scatter(x[:,0],t,c=z[:,0]+2,cmap=cmap,vmin=-2, vmax=3,s=10)
        tt = np.ones(ele_x[magma>melt_thod].shape)*frame*0.2
        if len(ele_x[magma>melt_thod])!=0:
            kk = np.ones(len(ele_x[magma>melt_thod]))*frame*0.2
            ax.scatter(ele_x[magma>melt_thod],kk,c='r',alpha = 0.2,s=5)
        if len(sxx[sxx>50])!=0:
            tip_sxx = sxx[sxx>50][-1]
            tip_z[frame] = ele_z[np.where(sxx==tip_sxx)]
            tip_x[frame] = ele_x[np.where(sxx==tip_sxx)]
        tttime[frame] = fl.time[frame]
        ### eclogite
        for xx in range(len(ele_x)):
            if True in (phase[xx,:]==phase_eclogite):
                for zz in range(len(ele_x[0])):
                    if (phase[xx,zz]==phase_eclogite):
                        x_change_ec[frame]=ele_x[xx,zz]
                        z_change_ec[frame]=-ele_z[xx,zz]
                        break
                break
        
        ### serpentinite
        for xx in range(len(ele_x)-1,0,-1):
            if True in (phase[xx,:]==phase_serpentinite):
                for zz in range(len(ele_x[0])-1,0,-1):
                    if (phase[xx,zz]==phase_serpentinite):
                        x_change_serp[frame]=ele_x[xx,zz]
                        z_change_serp[frame]=-ele_z[xx,zz]
                        break
                break
        
        ### chlorite
        for xx in range(len(ele_x)-1,0,-1):
            if True in (phase[xx,:]==phase_hydratedmantle):
                for zz in range(len(ele_x[0])-1,0,-1):
                    if (phase[xx,zz]==phase_hydratedmantle):
                        x_change_cho[frame]=ele_x[xx,zz]
                        z_change_cho[frame]=-ele_z[xx,zz]
                        break
                break
    #ax.scatter(tip_x,tttime,c=newcolors[9],s=10,label = 'bending')
    ax.scatter(x_change_ec,tttime,c=newcolors[1],s=10,label = 'basalt-eclogite')
    ax.scatter(x_change_serp,tttime,c=newcolors[0],s=10,label = 'serpentinite-peridotite')
    ax.scatter(x_change_cho,tttime,c=newcolors[3],s=10,label = 'chlorite-peridotite')
    ax.scatter(np.max(ele_x[magma>1.5e-3]),frame*0.2,c='r',label = 'melt')
    # name=model+'_flatslab_time_len.txt'
    # time,length,depth=np.loadtxt(savepath+name).T
    # name=model+'_flat_time_len.txt'
    # time_flat,flength,fdepth=np.loadtxt(savepath+name).T
    
    # ---------------------- plot setting --------------------------
    #ax.set_aspect('equal')
    for aa in [ax]:
        aa.set_xlim(300,xmax)
        aa.tick_params(labelsize=23)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
        
    #ymajor_ticks = np.linspace(300,0,num=7)
    #ax.set_yticks(ymajor_ticks)
    ax.set_ylim(10,45)
    #ax2.set_ylim(-0.1,0.2)
    #ax3.set_ylim(-2,4)
    ax.legend(fontsize=14)
    #ax3.tick_params(labelcolor="#4198B9")
    ax.set_ylabel('time (Myr)',fontsize=20)
    ax.set_xlabel('distance (km)',fontsize=20)
    ax_cbin = fig4.add_axes([0.93, 0.125, 0.02, 0.757])
    
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='vertical')
    cb.set_label(label='bathymetry (km)', size=25)
    cb.ax.yaxis.set_label_position('right')
    cb.ax.tick_params(labelsize=20)
if plot_sxx:
    cmap = plt.cm.get_cmap('seismic')
    fig5, (ax)= plt.subplots(1,1,figsize=(8,10))
    for frame in range(1,end,2):

        x,z = fl.read_mesh(frame)
        elex,elez = flac.elem_coord(x,z)
        sxx = fl.read_sxx(frame)*100
        xt = elex[:,0]
        t = np.zeros(xt.shape)
        t[:] = frame*0.2
        cb_plot=ax.scatter(elex[:,0],t,c=sxx[:,2],cmap=cmap,vmin=-200, vmax=200,s=35)
    
    # ---------------------- plot setting --------------------------
    #ax.set_aspect('equal')
    for aa in [ax]:
        aa.set_xlim(300,xmax)
        aa.tick_params(labelsize=23)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
        
    #ymajor_ticks = np.linspace(300,0,num=7)
    #ax.set_yticks(ymajor_ticks)
    ax.set_ylim(10,45)
    #ax2.set_ylim(-0.1,0.2)
    #ax3.set_ylim(-2,4)
    #ax3.tick_params(labelcolor="#4198B9")
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax_cbin = fig5.add_axes([0.93, 0.125, 0.02, 0.757])
    
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='vertical')
    cb.set_label(label='Bathymetry (km)', size=25)
    cb.ax.yaxis.set_label_position('right')
    cb.ax.tick_params(labelsize=20)
    
    
if plot_s1_profile:  
    skip = (slice(None, None, 5), slice(None, None, 3))
    fig, (ax)= plt.subplots(1,1,figsize=(17,16))
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- plotting -------------------------
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    ax.plot(ele_x[:,0],-ele_z[:,0],c = 'k',lw=3,linestyle='dashed')
    sxx = fl.read_sxx(frame)
    sxz = fl.read_sxz(frame)
    szz = fl.read_szz(frame)
    s1,s3,s2 = compute_s1(sxx, szz, sxz).T
    cbsxx = plt.cm.get_cmap('seismic_r')
    cbsxx = plt.cm.get_cmap('bwr_r')
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx*100,cmap=cbsxx,vmin=-500,vmax=500,shading='gouraud')
    
    ax.quiver(ele_x[skip],-ele_z[skip],s1.T[skip],s3.T[skip],headwidth=0,color = 'green')#,
              #angles='xy', scale_units='xy', scale=0.1,)
    ax.set_title('sxx',fontsize=25)
    cax = plt.axes([0.945, 0.365, 0.01, 0.271])
    cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
    cc1.ax.tick_params(labelsize=20)
    ax.contour(x,-z,temp,levels =[400,600,800],linewidths=3,colors = '#696969')
    cc1.set_label(label='$\sigma_{xx}$ (MPa)', size=25)
    cc1.ax.yaxis.set_label_position('left')
    ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    # ---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(300,0,num=7)
    ax.set_yticks(ymajor_ticks)
    ax.set_ylim(-zmin,-zmax)
    ax.set_ylim(230,-10)
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(300,950)
    #xmajor_ticks = np.linspace(250,1000,num=7)
    #ax.set_xticks(xmajor_ticks)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
    xx,zz,xt = np.loadtxt(savepath+'Nazca_aa06_40.0_final_slab.txt').T 
    xx=xx[zz<0]
    zz=zz[zz<0]
    zz = fd.moving_window_smooth(zz,3)
    ax.plot(xx+shift,-zz,color='k',lw=5)
       
    xxm,zzm,xtm = np.loadtxt(savepath+'Nazca_aa06_40_final_moho_slab.txt').T
    xxm=xxm[zzm<0]
    zzm=zzm[zzm<0]
    zzm = fd.moving_window_smooth(zzm,3)
    ax.plot(xxm+shift-5,-zzm+2,color='k',lw=5)
    
    
    ax.set_xlabel('distance (km)',fontsize=26)
    ax.set_ylabel('depth (km)',fontsize=26)
    
    ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
    fig.savefig(figpath+model+'frame_'+str(frame)+'_sxx_gourand.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))


if plot_sII:
    cmap = plt.cm.get_cmap('seismic')
    fig6, (ax)= plt.subplots(1,1,figsize=(8,10))
    for frame in range(1,end,2):

        x,z = fl.read_mesh(frame)
        elex,elez = flac.elem_coord(x,z)
        sII = fl.read_sII(frame)
        xt = elex[:,0]
        t = np.zeros(xt.shape)
        t[:] = frame*0.2
        cb_plot=ax.scatter(elex[:,0],t,c=sII[:,0],cmap=cmap,vmin=0, vmax=1,s=35)
    
    # ---------------------- plot setting --------------------------
    #ax.set_aspect('equal')
    for aa in [ax]:
        aa.set_xlim(300,xmax)
        aa.tick_params(labelsize=23)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
        
    #ymajor_ticks = np.linspace(300,0,num=7)
    #ax.set_yticks(ymajor_ticks)
    ax.set_ylim(10,45)
    #ax2.set_ylim(-0.1,0.2)
    #ax3.set_ylim(-2,4)
    #ax3.tick_params(labelcolor="#4198B9")
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax_cbin = fig6.add_axes([0.93, 0.125, 0.02, 0.757])
    
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='vertical')
    cb.set_label(label='Bathymetry (km)', size=25)
    cb.ax.yaxis.set_label_position('right')
    cb.ax.tick_params(labelsize=20)

if plot_strength_litho:
    fig7, (ax) = plt.subplots(1,1,figsize=(8,10))
    cmap = plt.cm.get_cmap('gist_earth')
    #os.chdir('/Users/chingchen/Desktop/pythonFLAC/')
    for frame in range(1,end,2):
        # ----------------------------- initial setup ---------------------------------
        phase = fl.read_phase(frame)
        x,z = fl.read_mesh(frame)
        elex,elez = flac.elem_coord(x,z)
        vis = 10**(fl.read_visc(frame))
        srII = 10**(fl.read_srII(frame))
        density = fl.read_density(frame)
        x_strength = np.zeros(len(elex[:,0]))
        for xx in range(0,len(elex[:,0])):
            #ind_continent = ((phase[xx,0]==phase_lowercrust)+(phase[xx,0]==phase_uppercrust))
            #if ind_continent :
                #print(xx)
            pressure=0
            total_strength = 0
            for zz in range(10,20):
                if phase[xx,zz]==phase_sediment or phase[xx,zz]==phase_serpentinite:
                    phi = np.tan(math.radians(3))
                    coh = 4e6
                else:
                    phi = np.tan(math.radians(30))
                    coh = 4e7
                layer = -elez[xx,zz]*1e3-(-elez[xx,zz-1]*1e3)
                if elez[0,zz] < -100:
                    break
                
                pp =  density[xx,zz]*9.8*layer
                #pressure =pp
                pressure +=pp
                strvis = 2 * vis[xx,zz]*srII[xx,zz]
                #print('-----------------',elez[xx,zz],'-----------------')
                strpls = coh + phi*pressure
                final_str = np.amin((strvis,strpls),axis=0) #Pa
                total_strength +=final_str
                #total_strength =final_str
                #print(final_str,total_strength,' Gpa')
            x_strength[xx] = total_strength
        xt = elex[:,0]
        t = np.zeros(xt.shape)
        t[:] = frame*0.2
        cb_plot=ax.scatter(elex[:,0],t,c=x_strength,cmap=cmap,s=35)
        
    # ---------------------- plot setting --------------------------
    #ax.set_aspect('equal')
    for aa in [ax]:
        aa.set_xlim(300,xmax)
        aa.tick_params(labelsize=23)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
        
    #ymajor_ticks = np.linspace(300,0,num=7)
    #ax.set_yticks(ymajor_ticks)
    ax.set_ylim(10,45)
    #ax2.set_ylim(-0.1,0.2)
    #ax3.set_ylim(-2,4)
    #ax3.tick_params(labelcolor="#4198B9")
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    ax_cbin = fig7.add_axes([0.93, 0.125, 0.02, 0.757])
    
    cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='vertical')
    cb.set_label(label='Bathymetry (km)', size=25)
    cb.ax.yaxis.set_label_position('right')
    cb.ax.tick_params(labelsize=20)