#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 19:20:27 2023

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

png             = 0
pdf             = 0

### plot
plot_topo_and_xvel=0
plot_topo_and_zvel=0
plot_s1_profile = 1
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


xmin,xmax=500,900
zmin,zmax=-120,10
model = 'Ref_Cocos'
shift = 550


frame = 51
frame = 11
frame = 190
frame = 126
frame = 201
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

    
g=10


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
#xmajor_ticks = np.linspace(250,1000,num=7)
#ax.set_xticks(xmajor_ticks)

ax.set_title('Time '+str(np.round(fl.time[frame-1],1))+' Myr',fontsize=20)
xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T

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
    cbsxx=ax.pcolormesh(ele_x,-ele_z,sxx*100,cmap=cbsxx,vmin=-2000,vmax=2000,shading='gouraud')
    
    ax.quiver(ele_x[skip],-ele_z[skip],s1.T[skip],s3.T[skip],headwidth=0)#,
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
    #fig.savefig(figpath+model+'frame_'+str(frame)+'_sxx_gourand.pdf')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
    # print("--- %s seconds ---" % (time.time() - start_time))
