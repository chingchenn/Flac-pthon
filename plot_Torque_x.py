#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 10:36:45 2022

@author: chingchen
"""

import numpy as np
import matplotlib.pyplot as plt
import function_for_flac as fd
plt.rcParams["font.family"] = "Times New Roman"
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
figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
gif = 0
mp4 = 0
end=150
model = 'Cocos_a0807'
# model = 'Nazca_a0706'
if model=='Nazca_a0706':
    xmin,xmax=0,200
    ymin,ymax=-5e18,5e18
    ymin3,ymax3=-1e17,1e17
    ymin4,ymax4=55,85
if model=='Nazca_a0702':
    xmin,xmax=0,800
    ymin,ymax=-5e18,3e19
    ymin3,ymax3=-1e17,5e17
    ymin4,ymax4=0,50
if model=='Ref_Cocos':
    xmin,xmax=0,300
    ymin,ymax=-1e18,5e18
    ymin3,ymax3=-1e17,1e17
    ymin4,ymax4=0,80
if model=='Cocos_a0807':
    xmin,xmax=0,150
    ymin,ymax=-1e18,5e18
    ymin3,ymax3=-1e17,1e17
    ymin4,ymax4=0,80
bwith = 3
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T

#for i in  range(5,end):
for i in  range(10,120):
    if i < 10:
        qq = '00'+str(i)
    elif i < 100 and i >=10:
        qq = '0'+str(i)
    else:
        qq=str(i)
    xg,fg, fgs = np.loadtxt(savepath+model+'_gravityx_'+str(qq)+'.txt').T
    xs,fs, fss,beta = np.loadtxt(savepath+model+'_suctionx_'+str(qq)+'.txt').T
    fig2, (ax2,ax3,ax4)= plt.subplots(3,1,figsize=(12,12))
    ax2.plot(xg-trench_x[i],fg ,c ='#c06c84',label='slab pull (N)',lw=3) 
    ax2.plot(xs-trench_x[i],fs,c="#FF9900",label='suction force (N)',lw=3)
    ax2.axhline(y=0, color='#849DAB', linestyle='--',lw=2)
    # print(np.min(-beta))
    # print(min(np.arccos(-beta)*180/np.pi))
    ax2.grid()
    # print(np.max(fs))
    fss = fd.moving_window_smooth(fss, 5)
    ax3.plot(xg-trench_x[i],fgs ,c ='#c06c84',label='slab pull (N)',lw=3) 
    ax3.plot(xs-trench_x[i],fss ,c ='#FF9900',label='suction force (N)',lw=3) 
    ax4.scatter(xs-trench_x[i],beta,c = 'darkgreen')
    ax2.set_ylim(ymin,ymax)
    ax3.set_ylim(ymin3,ymax3)
    ax4.set_ylim(ymin4,ymax4)
    ax2.set_title(model+' '+str(round(i*0.2,1))+' Myr',fontsize=24)
    ax3.set_xlabel('Time (Myr)',fontsize=26)
    ax2.set_ylabel('Torque (N)',fontsize=26)
    for axx in [ax2,ax3,ax4]:
        axx.set_xlim(xmin,xmax)
        axx.tick_params(axis='x', labelsize=26)
        axx.tick_params(axis='y', labelsize=26)
        axx.spines['bottom'].set_linewidth(bwith)
        axx.spines['top'].set_linewidth(bwith)
        axx.spines['right'].set_linewidth(bwith)
        axx.spines['left'].set_linewidth(bwith)
    
    ax2.legend(fontsize=16,loc='upper left')
    #fig2.savefig(figpath+model+'torquex'+str(qq)+'.png')
    #plt.close(fig2)

#-----------------------------creat GIF-----------------------------------------
if gif: 
    from PIL import Image
     
    # Create the frames
    frames = []
    for i in  range(5,end):
        if i < 10:
            qq = '00'+str(i)
        elif i < 100 and i >=10:
            qq = '0'+str(i)
        else:
            qq=str(i)
        img=figpath+model+'torquex'+str(qq)+'.png'
        new_frame = Image.open(img)
        frames.append(new_frame)
     
    # Save into a GIF file that loops forever
    frames[0].save(figpath+'png_to_gif.gif', format='GIF', append_images=frames[1:], 
                   save_all=True, duration=75, loop=0)
    
#-----------------------------creat mp4-----------------------------------------    
if mp4:
    import moviepy.editor as mp
    clip = mp.VideoFileClip(figpath+'png_to_gif.gif')
    # clip.write_videofile(figpath+'allfield_150'+model+".mp4")
    clip.write_videofile(figpath+'torque'+model+".mp4")



