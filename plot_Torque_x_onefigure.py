#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 21:02:00 2022

@author: chingchen
"""

import os
import flac
import numpy as np
import matplotlib
from matplotlib import cm
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
figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
gif = 0
mp4 = 0
end=150
model = 'Nazca_a0702'
# model = 'Nazca_a0706'
if model=='Nazca_a0706':
    xmin,xmax=0,300
    ymin,ymax=-5e18,5e18
    ymin3,ymax3=-1e17,1e17
    ymin4,ymax4=0,85
if model=='Nazca_a0702':
    xmin,xmax=200,1000
    ymin,ymax=-0.5,3.5
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



os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
cm = plt.cm.get_cmap('RdYlBu_r')

fig, (ax)= plt.subplots(2,2,figsize=(20,10))
ww = '126'
i = int(ww)
x, z = fl.read_mesh(i)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(i)
density = fl.read_density(i)
ax[0,0].pcolormesh(x,-z,phase,cmap=phase15,vmin=1, vmax=20)
ax[0,0].set_ylim(200,-5)
ax[0,0].set_aspect('equal')
temp = fl.read_temperature(i)
ax[0,0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
xg,fg, fgs = np.loadtxt(savepath+model+'_gravityx_'+str(ww)+'.txt').T
xs,fs, fss,beta = np.loadtxt(savepath+model+'_suctionx_'+str(ww)+'.txt').T
# ax[2,0].plot(xg-trench_x[i],fg ,c ='#c06c84',label='slab pull (N)',lw=3) 
# ax[2,0].plot(xs-trench_x[i],fs,c="#FF9900",label='suction force (N)',lw=3)
ax[1,0].plot(xg,fg/1e19 ,c ='#c06c84',label='slab pull (N)',lw=3) 
ax[1,0].plot(xs,fs/1e19,c="#FF9900",label='suction force (N)',lw=3)
ax[1,0].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,0].set_title(str(round(i/5,0))+' Myr',fontsize=28)

ww = '151'
i = int(ww)
x, z = fl.read_mesh(i)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(i)
density = fl.read_density(i)
ax[0,1].pcolormesh(x,-z,phase,cmap=phase15,vmin=1, vmax=20)
ax[0,1].set_ylim(200,-5)
ax[0,1].set_aspect('equal')
temp = fl.read_temperature(i)
ax[0,1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
xg,fg, fgs = np.loadtxt(savepath+model+'_gravityx_'+str(ww)+'.txt').T
xs,fs, fss,beta = np.loadtxt(savepath+model+'_suctionx_'+str(ww)+'.txt').T
# ax[3,0].plot(xg-trench_x[i],fg ,c ='#c06c84',label='slab pull (N)',lw=3) 
# ax[3,0].plot(xs-trench_x[i],fs,c="#FF9900",label='suction force (N)',lw=3)
ax[1,1].plot(xg,fg/1e19 ,c ='#c06c84',label='slab pull (N)',lw=3) 
ax[1,1].plot(xs,fs/1e19,c="#FF9900",label='suction force (N)',lw=3)
ax[1,1].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,1].set_title(str(round(i/5,0))+' Myr',fontsize=28)
ax[1,1].grid()
ax[1,0].grid()
ax[1,0].set_ylabel('Torque (10$^{19}$ N)',fontsize=30)
ax[1,0].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,1].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
for yy in range(len(ax)):
    for qq in range(len(ax[0])):  
        ax[1,qq].set_ylim(ymin,ymax)
        ax[yy,qq].set_xlim(xmin,xmax)
        # ax[yy,0].set_ylabel('Depth (km)',fontsize=20)
        ax[-1,qq].set_xlabel('Distance (km)',fontsize=30)
        bwith = 3
        ax[yy,qq].spines['bottom'].set_linewidth(bwith)
        ax[yy,qq].spines['top'].set_linewidth(bwith)
        ax[yy,qq].spines['right'].set_linewidth(bwith)
        ax[yy,qq].spines['left'].set_linewidth(bwith)
        ax[yy,qq].tick_params(axis='x', labelsize=30)
        ax[yy,qq].tick_params(axis='y', labelsize=30)
fig.tight_layout()
# fig.savefig(figpath+model+'frame_'+ww+'_Torquex.pdf')

fig2, (ax)= plt.subplots(2,2,figsize=(20,10))
ww = '176'
i = int(ww)
x, z = fl.read_mesh(i)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(i)
density = fl.read_density(i)
ax[0,0].pcolormesh(x,-z,phase,cmap=phase15,vmin=1, vmax=20)
ax[0,0].set_ylim(200,-5)
ax[0,0].set_aspect('equal')
temp = fl.read_temperature(i)
ax[0,0].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
xg,fg, fgs = np.loadtxt(savepath+model+'_gravityx_'+str(ww)+'.txt').T
xs,fs, fss,beta = np.loadtxt(savepath+model+'_suctionx_'+str(ww)+'.txt').T
# ax[2,0].plot(xg-trench_x[i],fg ,c ='#c06c84',label='slab pull (N)',lw=3) 
# ax[2,0].plot(xs-trench_x[i],fs,c="#FF9900",label='suction force (N)',lw=3)
ax[1,0].plot(xg,fg/1e19 ,c ='#c06c84',label='slab pull (N)',lw=3) 
ax[1,0].plot(xs,fs/1e19,c="#FF9900",label='suction force (N)',lw=3)
ax[1,0].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,0].set_title(str(round(i/5,0))+' Myr',fontsize=28)

ww = '201'
i = int(ww)
x, z = fl.read_mesh(i)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(i)
density = fl.read_density(i)
ax[0,1].pcolormesh(x,-z,phase,cmap=phase15,vmin=1, vmax=20)
ax[0,1].set_ylim(200,-5)
ax[0,1].set_aspect('equal')
temp = fl.read_temperature(i)
ax[0,1].contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
xg,fg, fgs = np.loadtxt(savepath+model+'_gravityx_'+str(ww)+'.txt').T
xs,fs, fss,beta = np.loadtxt(savepath+model+'_suctionx_'+str(ww)+'.txt').T
# ax[3,0].plot(xg-trench_x[i],fg ,c ='#c06c84',label='slab pull (N)',lw=3) 
# ax[3,0].plot(xs-trench_x[i],fs,c="#FF9900",label='suction force (N)',lw=3)
ax[1,1].plot(xg,fg/1e19 ,c ='#c06c84',label='slab pull (N)',lw=3) 
ax[1,1].plot(xs,fs/1e19,c="#FF9900",label='suction force (N)',lw=3)
ax[1,1].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,1].set_title(str(round(i/5,0))+' Myr',fontsize=28)
ax[1,1].grid()
ax[1,0].grid()
ax[1,0].set_ylabel('Torque (10$^{19}$ N)',fontsize=30)
ax[1,0].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
ax[1,1].axhline(y=0, color='#849DAB', linestyle='--',lw=2)
for yy in range(len(ax)):
    for qq in range(len(ax[0])):
        ax[1,qq].set_ylim(ymin,ymax)
        ax[yy,qq].set_xlim(xmin,xmax)
        # ax[yy,0].set_ylabel('Depth (km)',fontsize=20)
        ax[-1,qq].set_xlabel('Distance (km)',fontsize=30)
        bwith = 3
        ax[yy,qq].spines['bottom'].set_linewidth(bwith)
        ax[yy,qq].spines['top'].set_linewidth(bwith)
        ax[yy,qq].spines['right'].set_linewidth(bwith)
        ax[yy,qq].spines['left'].set_linewidth(bwith)
        # ax[yy,0].spines['top'].set_visible(False)
        ax[yy,qq].tick_params(axis='x', labelsize=30)
        ax[yy,qq].tick_params(axis='y', labelsize=30)
fig2.tight_layout()
# fig.savefig(figpath+model+'frame_'+ww+'_Torquex.pdf')