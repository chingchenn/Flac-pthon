#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 09:39:22 2022

@author: ji-chingchen
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


plt.rcParams["font.family"] = "Times New Roman"
path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Volumes/SSD500/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'

fig, axes = plt.subplots(1, 1, figsize=(15, 1))
fig.subplots_adjust(wspace=2.5)
labelsize=40

cmap1 = copy.copy(cm.jet)
norm1 = mcolors.Normalize(vmin=20, vmax=27)
im1 = cm.ScalarMappable(norm=norm1, cmap=cmap1)
cbar1 = fig.colorbar(im1, cax=axes,ticks=np.linspace(20, 27, 8),orientation='horizontal')
cbar1.ax.tick_params(labelsize=labelsize) 
# cbar1.set_label('# of contacts', rotation=270,fontsize = 20)
fig2, axes2 = plt.subplots(1,1, figsize=(15, 1))
bins = [0,200,400,600,800,1000,1200]
nbin = len(bins) - 1
cmap4 = cm.get_cmap('rainbow', nbin)
norm4 = mcolors.BoundaryNorm(bins, nbin)
im4 = cm.ScalarMappable(norm=norm4, cmap=cmap4)
cbar4 = fig.colorbar(im4, cax=axes2,orientation='horizontal')
cbar4.ax.tick_params(labelsize=labelsize) 
fig.savefig(figpath+'viscosity.pdf')
fig2.savefig(figpath+'temperature.pdf')
