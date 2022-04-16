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
#matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
from matplotlib import animation
import Main_snapshot as Ms
plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'

# model = sys.argv[1]
model = 'Chi02'
# frame = int(sys.argv[2])
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
time=fl.time
gif = 0
mp4 = 0

#------------------------------------------------------------------------------
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#F87431","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)

for i in range(95,96):
    fig, (ax,ax2)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
    x,z,ele_x,ele_z,phase,temp,ztop=Ms.plot_snapshot(i)
    ax.scatter(ele_x,-ele_z,c = phase,cmap = phase19,vmax=19,vmin=1,s=150)
    ax.contour(x,-z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
    
    #--------------------------------------------------------------------------
    ele_x,ele_z,vis,ztop = Ms.get_vis(i)
    cc = plt.cm.get_cmap('rainbow')
    cb_plot=ax2.scatter(ele_x,-ele_z,c=vis,cmap=cc,vmin=20, vmax=27,s=150)
    # ax_cbin = fig.add_axes([0.63, 0.35, 0.23, 0.03])
    # cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    
    # ax_cbin.set_title('Pa s',fontsize=20)
    
    
    # ax2.set_title('Viscous '+str(model)+' at '+str(round(fl.time[i-1],1))+' Myr',fontsize=24)

    

    # ax.set_ylabel('Depth (km)',fontsize=20)
    # ax.set_xlabel('Distance (km)',fontsize=20)
    # ax2.set_ylabel('Depth (km)',fontsize=20)
    # ax2.set_xlabel('Distance (km)',fontsize=20)
    ax.set_aspect('equal')
    ax2.set_aspect('equal')
    
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.tick_params(axis='x', labelsize=23)
    ax.tick_params(axis='y', labelsize=23)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax2.tick_params(axis='x', labelsize=23)
    ax2.tick_params(axis='y', labelsize=23)
    ymajor_ticks = np.linspace(200,0,num=5)
    ax.set_yticks(ymajor_ticks)
    ax2.set_yticks(ymajor_ticks)
    xmajor_ticks = np.linspace(250,850,num=5)
    ax.set_xticks(xmajor_ticks)
    ax2.set_xticks(xmajor_ticks)
    ax.set_xlim(250,850)
    ax.set_ylim(200,-30)
    ax.set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
    ax2.set_xlim(250,850)
    ax2.set_ylim(200,-30)
    if i < 10:
        qq = '00'+str(i)
    elif i < 100 and i >=10:
        qq = '0'+str(i)
    else:
        qq=str(i)
    plt.show()
    fig.savefig(path+model+'/frame_'+qq+'_phase_vis.png')
    fig.gca()
    plt.close(fig)

#-----------------------------creat GIF-----------------------------------------
if gif: 
    from PIL import Image
    import glob
     
    # Create the frames
    frames = []
    imgs = glob.glob(path+model+"/*_phase_vis.png")
    for i in imgs:
        new_frame = Image.open(i)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save('png_to_gif.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=300, loop=0)
    
#-----------------------------creat mp4-----------------------------------------    
if mp4:
    import moviepy.editor as mp
    clip = mp.VideoFileClip("png_to_gif.gif")
    clip.write_videofile(model+".mp4")

