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

model=sys.argv[1]
name='melting_'+model
fig, (ax) = plt.subplots(1,1,figsize=(18,9))
time,phase_p4,phase_p9,phase_p10 = np.loadtxt('/home/jiching/geoflac/data/'+name+'.txt').T
ax.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
ax.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments')
ax.set_xlim(0,30)
#ax.grid()
ax.tick_params(axis='x', labelsize=26)
ax.tick_params(axis='y', labelsize=26)
ax.legend(fontsize=25)
#ax.set_title('Model : '+model,fontsize=25)
ax.set_xlabel('Time (Myr)',fontsize=26)
ax.set_ylabel('molten rocks (km3/km)',fontsize=26)
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
figpath = '/home/jiching/geoflac/figure/'
fig.savefig(figpath+model+'_single_bar_plot_melting.png')
