#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 11:04:51 2022

@author: ji-chingchen
"""

import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
import Main_snapshot as Ms
# plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'

model = sys.argv[1]
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
time=fl.time
plotting_png = 1
plotting_vx = 0
plotting_vz = 0
mp4 = 1
labelsize = 26
if not os.path.isdir(path+model+'/frame_plot'):
    os.mkdir(path+model+'/frame_plot')
#------------------------------------------------------------------------------
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)
if plotting_png:
    for i in range(1,end,10):
        fig, (ax,ax2)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
        x,z,ele_x,ele_z,phase,temp,ztop=Ms.plot_snapshot(i)
        xm, zm, age, ph, idd, a1, a2, ntriag = fl.read_markers(i)
        #ax.scatter(ele_x,-ele_z,c = phase,cmap = phase19,vmax=19,vmin=1,s=150)
        ax.scatter(xm,-zm,c = ph,cmap = phase19,vmax=19,vmin=1,s=5)
        ax.contour(x,-z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
    #--------------------------------------------------------------------------
        ele_x,ele_z,vis,ztop = Ms.get_vis(i)
        cc = plt.cm.get_cmap('jet')
        cb_plot=ax2.scatter(ele_x,-ele_z,c=vis,cmap=cc,vmin=20, vmax=27,s=150)
        ax.set_aspect('equal')
        ax2.set_aspect('equal')
     
        bwith = 3
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.tick_params(axis='x', labelsize=labelsize)
        ax.tick_params(axis='y', labelsize=labelsize)
        ax2.spines['bottom'].set_linewidth(bwith)
        ax2.spines['top'].set_linewidth(bwith)
        ax2.spines['right'].set_linewidth(bwith)
        ax2.spines['left'].set_linewidth(bwith)
        ax2.tick_params(axis='x', labelsize=labelsize)
        ax2.tick_params(axis='y', labelsize=labelsize)
        ymajor_ticks = np.linspace(200,0,num=7)
        ax.set_yticks(ymajor_ticks)
        ax2.set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1200,num=6)
        ax.set_xticks(xmajor_ticks)
        ax2.set_xticks(xmajor_ticks)
        ax.set_xlim(250,1200)
        ax.set_ylim(200,-30)
        ax.set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
        ax2.set_xlim(250,1200)
        ax2.set_ylim(250,-30)
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        fig.savefig(path+model+'/frame_plot/frame_'+qq+'_phase_vis.png')
        fig.gca()
        plt.close(fig)
        print('----- finish figure '+qq+' -----')
if plotting_vx:
    for i in range(1,end):
        fig,(ax)=plt.subplots(1,1,figsize=(12,8))
        cc = plt.cm.get_cmap('BrBG')
        x,z, = fl.read_mesh(i)
        x,z,ele_x,ele_z,phase,temp,ztop=Ms.plot_snapshot(i)
        vx,vz = fl.read_vel(i) # cm/yr
        ref_vx = vx[370,39]
        cb_plot = ax.scatter(x,z,c=(vx-ref_vx),cmap = cc,vmin = -3, vmax = 3)
        ax_cbin = fig.add_axes([0.2,0.4,0.15,0.01])
        cb = fig.colorbar(cb_plot,cax = ax_cbin,orientation = 'horizontal')
        ax.contour(x,z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        ax.set_ylabel('Depth (km)',fontsize=20)
        ax.set_xlabel('Distance (km)',fontsize=20)
        ax_cbin.set_title('cm/yr',fontsize=20)
        ax.set_aspect('equal')
        ax.set_xlim(0,1200)
        ax.set_ylim(-300,20)
        bwith = 3
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        ax.set_title('horizontal velocity (cm/yr) '+str(model)+' at '+str(round(fl.time[i-1],1))+' Myr',fontsize=24)
        fig.savefig(path+model+'/frame_plot/frame_'+qq+'_Vx.png')
        print('----- finish figure '+qq+' -----')
if plotting_vz:
    for i in range(1,end):
        fig,(ax)=plt.subplots(1,1,figsize=(12,8))
        cc = plt.cm.get_cmap('BrBG')
        x,z, = fl.read_mesh(i)
        x,z,ele_x,ele_z,phase,temp,ztop=Ms.plot_snapshot(i)
        vx,vz = fl.read_vel(i) # cm/yr
        ref_vz = vz[370,39]
        cb_plot = ax.scatter(x,z,c=(vz-ref_vz),cmap = cc,vmin = -5, vmax = 5)
        ax_cbin = fig.add_axes([0.2,0.4,0.15,0.01])
        cb = fig.colorbar(cb_plot,cax = ax_cbin,orientation = 'horizontal')
        ax.contour(x,z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        ax.set_ylabel('Depth (km)',fontsize=20)
        ax.set_xlabel('Distance (km)',fontsize=20)
        ax_cbin.set_title('cm/yr',fontsize=20)
        ax.set_aspect('equal')
        ax.set_xlim(0,1200)
        ax.set_ylim(-300,20)
        bwith = 3
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        ax.set_title('vertical velocity (cm/yr) '+str(model)+' at '+str(round(fl.time[i-1],1))+' Myr',fontsize=24)
        fig.savefig(path+model+'/frame_plot/frame_'+qq+'_Vz.png')
        print('----- finish figure '+qq+' -----')
#-----------------------------creat GIF-----------------------------------------
if mp4: 
    from PIL import Image
    import glob
     
    # Create the frames
    frames = []
    for i in  range(1,end,10):
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        if plotting_png:
            name = '_phase_vis'
        elif plotting_vx:
            name = '_Vx'
        elif plotting_vz:
            name = '_Vz'
        else:
            name = '_phase_vis'
        img=path+model+'/frame_plot/frame_'+qq+name+'.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(path+model+'/frame_plot/'+name+'.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=80, loop=0)
#-----------------------------creat mp4-----------------------------------------    
    import moviepy.editor as mp
    clip = mp.VideoFileClip(path+model+"/frame_plot/"+name+".gif")
    #clip.write_videofile(figpath+'phase_vis_'+model+".mp4")
    clip.write_videofile(figpath+model+name+".mp4")
