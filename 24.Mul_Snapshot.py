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
# matplotlib.use('Agg')
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
path = '/Users/chingchen/Desktop/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
model_list = ['Nazca_a0702','Nazca_a0706']
# model_list = ['Nazca_a0634','Nazca_a0636']
#model_list=['Ref_Cocos','Cocos_a0807']
plotting_vis = 1
plotting_pre = 0
gif = 1
mp4 = 1
end=150
model = 'MUL'
chamber_limit = 5e-3
###------------------------------------------------------------------------------
def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre
#------------------------------------------------------------------------------
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)
for i in range(1,end+1):
#for i in range(end,end+1):
    if plotting_vis ==0:
        break
    fig, (ax)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
    os.chdir(path+model_list[0])
    fl = flac.Flac()
    ax[0].set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
    for kk,model in enumerate(model_list):
        os.chdir(path+model)
        fl = flac.Flac()
        x,z = fl.read_mesh(i)
        temp = fl.read_temperature(i)
        
        magma_chamber = fl.read_fmagma(i)*100
        ax[kk].contour(x,-z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        x,z,ele_x,ele_z,vis,ztop = Ms.get_vis(i)
        cc = plt.cm.get_cmap('jet')
        cb_plot=ax[kk].pcolormesh(x,-z,vis,cmap = cc,vmax=24,vmin=20)
        ax[kk].scatter(ele_x[magma_chamber > chamber_limit],-ele_z[magma_chamber > chamber_limit],magma_chamber[magma_chamber > chamber_limit]*1e3,c = 'w',alpha = 0.2)
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
    # fig.savefig(figpath+'frame_'+qq+'_compare_vis.png')
    fig.gca()
    # plt.close(fig)
for i in range(1,end+1,20):
    if plotting_pre ==0:
        break
    fig, (ax)= plt.subplots(2,1,figsize=(20,16),clear = True,gridspec_kw={'height_ratios':[1,1]})
    os.chdir(path+model_list[0])
    fl = flac.Flac()
    ax[0].set_title(str(round(fl.time[i-1],1))+' Myr',fontsize=36)
    for kk,model in enumerate(model_list):
        os.chdir(path+model)
        fl = flac.Flac()
        x,z = fl.read_mesh(i)
        temp = fl.read_temperature(i)
        x,z,new_pre = dynamics_pressure(i)
        ck = plt.cm.get_cmap('RdYlBu_r')
        cbpre=ax[kk].pcolormesh(x,-z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200)
        cax = plt.axes([0.945, 0.385, 0.01, 0.231])
        cc3=fig.colorbar(cbpre, ax=ax,cax=cax)
        cc3.set_label(label='Pressure (MPa)', size=25)
        cc3.ax.tick_params(labelsize=20)
        cc3.ax.yaxis.set_label_position('left')
        ax[kk].contour(x,-z,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
        
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
    fig.savefig(figpath+'frame_'+qq+'_compare_pre.png')
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
        if plotting_pre:
            img=figpath+'frame_'+qq+'_compare_pre.png'
        if plotting_vis:
            img=figpath+'frame_'+qq+'_compare_vis.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(figpath+'frame_'+qq+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=60, loop=0)
    
#-----------------------------creat mp4-----------------------------------------    
if mp4:
    import moviepy.editor as mp
    clip = mp.VideoFileClip(figpath+'frame_'+qq+'png_to_gif.gif')
    if plotting_pre:
        clip.write_videofile(figpath+'pre_'+model+".mp4")
    if plotting_vis:
        clip.write_videofile(figpath+'vis_'+model+".mp4")
    
