#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 16:15:30 2023

@author: chingchen
"""

import flac
import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Helvetica"
#=========================setting=============================
path = '/Users/chingchen/Desktop/model/'
savepath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation/'

bwith=3
fontsize = 30

fig6, (ax1)= plt.subplots(1,1,figsize=(15,7))  
model='Nazca_aa06'
model = 'Nazca_v2_07'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
#--------------------------------- FIG1 melting phase -------------------------
name='melting_'+model
time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
total = (phase_p10+phase_p3+phase_p4)
ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#4169E1',label='eclogite')
#------------------------------figure setting---------------------------------
name='Nazca_a0702'+'_flatslab_time_len.txt'
time_flat,length,depth=np.loadtxt(savepath+name).T
for aaa in [ax1]:
    aaa.tick_params(labelsize=fontsize)
    aaa.grid()
    aaa.set_xlim(10,40)
    aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
    for axis in ['top','bottom','left','right']:
        aaa.spines[axis].set_linewidth(bwith)
ax1.set_ylim(0,6)    
ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
ax1.set_xlabel('time (Myr)',fontsize=fontsize)
fig6.tight_layout()
ax1.legend(fontsize=fontsize,facecolor='white')
#fig6.savefig(figpath+'fig3_v3.pdf')