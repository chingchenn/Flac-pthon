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
plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
path = '/scratch2/jiching/04model/'
#path = '/Users/chingchen/Desktop/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/scratch2/jiching/data/'
#savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
#figpath = '/Users/chingchen/Desktop/figure/'
# model_list = ['Nazca_a0702','Nazca_a0706']
model_list = ['Nazca_a0634','Nazca_a0636']
plotting_png = 1
gif = 0
mp4 = 0
end=150
model = 'MUL'
chamber_limit = 5e-3
#------------------------------------------------------------------------------
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)
for i in range(1,end+1,10):
# for i in range(end,end+1):
    if plotting_png ==0:
        break
    fig, (ax)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
    os.chdir(path+model_list[0])
    fl = flac.Flac();end = fl.nrec
    ax[0].set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
    for kk,model in enumerate(model_list):
        os.chdir(path+model)
        fl = flac.Flac();end = fl.nrec
        x,z = fl.read_mesh(i)
        temp = fl.read_temperature(i)
        magma_chamber = fl.read_fmagma(i)*100
        print(np.max(magma_chamber))
        ax[kk].contour(x,-z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        x,z,ele_x,ele_z,vis,ztop = Ms.get_vis(i)
        cc = plt.cm.get_cmap('jet')
        cb_plot=ax[kk].pcolormesh(x,-z,vis,cmap = cc,vmax=27,vmin=20)
        ax[kk].scatter(ele_x[magma_chamber > chamber_limit],-ele_z[magma_chamber > chamber_limit],magma_chamber[magma_chamber > chamber_limit]*1e5,c = 'w')
        ax[kk].set_aspect('equal')
        bwith = 3
        ax[kk].spines['bottom'].set_linewidth(bwith)
        ax[kk].spines['top'].set_linewidth(bwith)
        ax[kk].spines['right'].set_linewidth(bwith)
        ax[kk].spines['left'].set_linewidth(bwith)
        ax[kk].tick_params(axis='x', labelsize=23)
        ax[kk].tick_params(axis='y', labelsize=23)
        ymajor_ticks = np.linspace(200,0,num=5)
        ax[kk].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=6)
        ax[kk].set_xticks(xmajor_ticks)
        ax[kk].set_xlim(250,1000)
        ax[kk].set_ylim(200,-30)
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
    fig.savefig(savepath+'frame_'+qq+'_compare_vis.png')
    fig.gca()
    plt.close(fig)

#-----------------------------creat GIF-----------------------------------------
if gif: 
    from PIL import Image
    import glob
     
    # Create the frames
    frames = []
    for i in  range(1,end+1):
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        img=savepath+'frame_'+qq+'_compare_vis.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(savepath+'frame_'+qq+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=75, loop=0)
    
#-----------------------------creat mp4-----------------------------------------    
if mp4:
    import moviepy.editor as mp
    clip = mp.VideoFileClip(savepath+'frame_'+qq+'png_to_gif.gif')
    clip.write_videofile(figpath+'vis_'+model+".mp4")
    
