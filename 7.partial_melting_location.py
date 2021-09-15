#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:19:28 2021

@author: ji-chingchen
"""
import os
import flac
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt


model='w1022'
    # path = '/scratch2/jiching/'+model+'/'
    # path = '/home/jiching/geoflac/'+model+'/'
path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
os.chdir(path)
# fl = flac.Flac();end = fl.nrec
fl = flac.FlacFromVTK();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
print(model)    
    
def melting_phase():
    melt_num = np.zeros(end)
    phase_p4=np.zeros(end)
    phase_p9=np.zeros(end)
    phase_p10=np.zeros(end)
    po=np.zeros(end)
    for i in range(1,end):
        c=0;p9=0;p4=0;p10=0
        x, z = fl.read_mesh(i)
        mm=fl.read_fmelt(i)
        phase=fl.read_phase(i)
        for xx in range(len(mm)):
            for zz in range(len(mm[0])):
                if mm[xx,zz] != 0:
                    print()
                    if phase[xx,zz]==9:
                        p9 += 1
                    elif phase[xx,zz]==4:
                        p4 += 1
                    elif phase[xx,zz]==10:
                        p10 += 1
                    c +=1
        pk=c-p4-p9-p10
        melt_num[i]=c
        phase_p4[i]=p4
        phase_p9[i]=p9
        phase_p10[i]=p10
        po[i]=pk
    return phase_p4,phase_p9,phase_p10,po
time =fl.time
phase_p4,phase_p9,phase_p10,po=melting_phase()

fig, (ax) = plt.subplots(1,1,figsize=(18,12))
others=po
ax.bar(time,phase_p9,width=0.17,color='orange',label='serpentinite ')
ax.bar(time,phase_p4,bottom=phase_p9,width=0.17,color='seagreen',label='olivine')
ax.bar(time,phase_p10,bottom=phase_p9+phase_p4,width=0.17,color='tomato',label='sediments')
ax.bar(time,others,bottom=phase_p9+phase_p4+phase_p10,width=0.17,color='k',label='others')
ax.set_xlim(0,24)
ax.grid()
ax.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16 )
ax.set_title('Model : '+model,fontsize=25)
ax.set_xlabel('Time (Myr)',fontsize=20)
ax.set_ylabel('number of elements in phases',fontsize=20)
ax.legend(fontsize=25)
# fig.savefig(figpath+model+'_bar_plot_melting.png')


plt.hist([1, 2, 1], bins=[0, 1, 2, 3])

