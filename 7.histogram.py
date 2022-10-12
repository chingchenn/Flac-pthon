#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:19:28 2021

@author: ji-chingchen
"""
import os,sys
import flac
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

# model=sys.argv[1]
model = 'Nazca_0502'
#name='melting_'+model
#fig, (ax) = plt.subplots(1,1,figsize=(18,12))
#time,phase_p4,phase_p5,phase_p9,phase_p10,others = np.loadtxt('/home/jiching/geoflac/data/'+name+'.txt').T
#ax.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
#ax.bar(time,phase_p10+phase_p5,bottom=phase_p9+phase_p4,width=0.17,color='tomato',label='sediments')
#ax.bar(time,others,bottom=phase_p9+phase_p4+phase_p10+phase_p5,width=0.17,color='k',label='others')
#ax.set_xlim(0,30)
#ax.grid()
#ax.tick_params(axis='x', labelsize=16 )
#ax.tick_params(axis='y', labelsize=16 )
#ax.legend(fontsize=25)
#ax.set_title('Model : '+model,fontsize=25)
#ax.set_xlabel('Time (Myr)',fontsize=20)
#ax.set_ylabel('molten rocks (km3/km)',fontsize=20)
#bwith = 3
#ax.spines['bottom'].set_linewidth(bwith)
#ax.spines['top'].set_linewidth(bwith)
#ax.spines['right'].set_linewidth(bwith)
#ax.spines['left'].set_linewidth(bwith)
figpath = '/home/jiching/geoflac/figure/'
#fig.savefig(figpath+model+'_single_bar_plot_melting.png')


#plt.hist([1, 2, 1], bins=[0, 1, 2, 3])
path = '/Users/chingchen/Desktop/model/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
i=39
melt_num = np.zeros(end)
phase_p4=np.zeros(end)
phase_p5=np.zeros(end)
phase_p9=np.zeros(end)
phase_p10=np.zeros(end)
po=np.zeros(end)
c=0;p9=0;p4=0;p10=0;p5=0
x, z = fl.read_mesh(i)
mm=fl.read_fmelt(i)
phase=fl.read_phase(i)
area = fl.read_area(i)
pp=[]
for xx in range(len(mm)):
    for zz in range(len(mm[0])):
        if mm[xx,zz] != 0:
            print(mm)
            pp.append(phase[xx,zz]) 
            print(phase[xx,zz])
            if phase[xx,zz]==9:
                p9 += area[xx,zz]*mm[xx,zz]/1e6
            elif phase[xx,zz]==4:
                p4 +=area[xx,zz]*mm[xx,zz]/1e6
            elif phase[xx,zz]==10:
                p10 += area[xx,zz]*mm[xx,zz]/1e6
            elif phase[xx,zz]==5:
                p5 += area[xx,zz]*mm[xx,zz]/1e6
            c +=1

fig, (ax) = plt.subplots(1,1,figsize=(18,12))
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
num_bin=20
ax.hist(pp,num_bin)
# fig.savefig(figpath+model+str(i)+'_single_bar.png')
