#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 15:59:37 2022

@author: chingchen
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

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
path = '/Users/chingchen/Desktop/model/'

savepath='/home/jiching/geoflac/data/'
#savepath = '/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:\\OneDrive - 國立台灣大學/resarch/data/'
#savepath='D:/model/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
model_list=['Ref03','h0401','h0402','h0403']
model_list=['b0607k','b0608k','b0609k','b0611k','b0614k','b0615k']
model_list=['Cocos_a0646','Cocos_a0807']
#model_list=['Nazca_a0701','Nazca_a0702','Ref_Nazca','Nazca_a0704','Nazca_a0705']
names=['120km','130km','140km']
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
plt.rcParams["font.family"] = "Times New Roman"
bwith = 3





print('--- start plot geometry ---')
fig2, (ax2) = plt.subplots(1,2,figsize=(20,9))
model = 'Cocos_a0646'
# for kk,model in enumerate(model_list):
os.chdir(path+model)
fl = flac.Flac()
frame = 150
temp = fl.read_temperature(frame)
x,z, = fl.read_mesh(frame)
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
cx=ax2[0].contour(x-trench_x[-1],-z,temp,cmap = 'rainbow',levels =[550],linewidths=2.6)
xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
print(savepath+str(model)+'_final_slab.txt')
xx= fd.moving_window_smooth(xmean[xmean>0], 5)
ztop = fd.moving_window_smooth(ztop[xmean>0], 5)
xmean=xx
ax2[0].plot(xmean,-ztop,c='k',label=model,lw=3)
ax2[0].set_ylabel("Depth (km)",fontsize=28)


model = 'Cocos_a0807'
# for kk,model in enumerate(model_list):
os.chdir(path+model)
fl = flac.Flac()
frame = 150
temp = fl.read_temperature(frame)
x,z, = fl.read_mesh(frame)
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
# x,z,ele_x,ele_z,phase,temp,ztop = plot_snapshot(frame)
cx=ax2[1].contour(x-trench_x[-1],-z,temp,cmap = 'rainbow',levels =[550],linewidths=2.6)
xmean,ztop=np.loadtxt(savepath+str(model)+'_final_slab.txt').T
print(savepath+str(model)+'_final_slab.txt')
xx= fd.moving_window_smooth(xmean[xmean>0], 5)
ztop = fd.moving_window_smooth(ztop[xmean>0], 5)
xmean=xx
ax2[1].plot(xmean,-ztop,c='k',label=model,lw=3)
ax2[1].set_ylabel("Depth (km)",fontsize=28)

for qq in ax2:
    
    qq.set_xlabel("Distance relative to trench (km)",fontsize=28)
    # ax2.legend(fontsize=16)
    xmajor_ticks = np.linspace(0,100,num=5)
    qq.set_yticks(xmajor_ticks)
    qq.set_ylim(100,0)
    qq.set_xlim(0,300)
    qq.spines['bottom'].set_linewidth(bwith)
    qq.spines['top'].set_linewidth(bwith)
    qq.spines['right'].set_linewidth(bwith)
    qq.spines['left'].set_linewidth(bwith)
    qq.set_aspect('equal')
    qq.tick_params(axis='x', labelsize=26)
    qq.tick_params(axis='y', labelsize=26)
    qq.grid()
#fig2.savefig('D:\\OneDrive - 國立台灣大學/master03/Seminar/'+'multi_slab_analysis_'+model_list[0]+'_'+model_list[-1]+'.pdf')
fig2.savefig(figpath+'Curie_Point.pdf')
print('=========== DONE =============')