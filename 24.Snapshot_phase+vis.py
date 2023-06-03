#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 12:55:30 2022

@author: chingchen
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
import Main_snapshot as Ms
plt.rcParams["font.family"] = "Times New Roman"
import time
start_time = time.time()
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
path = '/Users/chingchen/Desktop/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
model_list = ['Nazca_a0516']
plotting_png = 1
gif = 0
mp4 = 0
end=150


# plot_frame=[31,61,91,121]
#------------------------------------------------------------------------------
# 

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)
cc = plt.cm.get_cmap('jet')
for model in model_list:
    for i in range(1,end+1,20):
    # for i in plot_frame:
        if plotting_png ==0:
            break
        fig, (ax)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
        os.chdir(path+model)
        fl = flac.Flac();end = fl.nrec
        ax[0].set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
        # for kk,model in enumerate(model_list):
            
        x,z = fl.read_mesh(i)
        temp = fl.read_temperature(i)
        ax[0].contour(x,-z,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=5)
        ax[1].contour(x,-z,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=5)
        x,z,ele_x,ele_z,phase,temp,ztop = Ms.plot_snapshot(i)
        ax[0].pcolormesh(x,-z,phase,cmap = phase19,vmax=20,vmin=1,shading = 'auto')
        # xx,zz,ph = Ms.inter_phase(i)
        print('finish interp frame'+str(i))
        # ax[0].pcolormesh(xx,-zz,ph,cmap = phase19,vmax=20,vmin=1,shading = 'auto')
        # ax[0].scatter(xx,-zz,c = ph,cmap = phase19,vmax=19,vmin=1,s=10)
        x,z,ele_x,ele_z,vis,ztop = Ms.get_vis(i)
        cb_plot=ax[1].pcolormesh(x,-z,vis,cmap = cc,vmax=27,vmin=20)
        for kk in range(len(ax)):
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
        # fig.savefig(figpath+model+'_'+'frame_'+qq+'_phase+vis_test.png')
        fig.gca()
        # plt.close(fig)

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
            img=figpath+model+'_'+'frame_'+qq+'_phase+vis.png'
            new_frame = Image.open(img)
            frames.append(new_frame)
         
        # Save into a GIF file that loops forever
        frames[0].save(savepath+model+'_'+'frame_'+qq+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                       save_all=True, duration=75, loop=0)
        
    #-----------------------------creat mp4-----------------------------------------    
    if mp4:
        import moviepy.editor as mp
        clip = mp.VideoFileClip(savepath+model+'_'+'frame_'+qq+'png_to_gif.gif')
        clip.write_videofile(figpath+'vis_'+model+".mp4")
    
print("--- %s seconds ---" % (time.time() - start_time))