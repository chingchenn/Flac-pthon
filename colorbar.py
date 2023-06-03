#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 23:46:15 2022

@author: ji-chingchen
"""

# import pylab as pl
# import numpy as np

# a = np.array([[0,1]])

# img = pl.imshow(a, cmap="jet")
# pl.gca().set_visible(False)
# cax = pl.axes([1, 0.2, 0.8, 0.6])
# pl.colorbar(orientation="horizontal", cax=cax)
# # pl.savefig("colorbar.pdf")


import matplotlib.pyplot as plt
import matplotlib as mpl
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/colorbar/'
# plt.figure(figsize=(9, 1.5))
fig, ax = plt.subplots(figsize=(2, 17))
# fig.subplots_adjust(bottom=0.5)

cmap = mpl.cm.rainbow

bounds =[200,400,600,800,1000,1200]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb2 = mpl.colorbar.ColorbarBase(ax,cmap=cmap,
                                 norm=norm,
                                orientation='vertical')
cb2.ax.tick_params(labelsize=40)
# fig.savefig(figpath+'temperature_ver.pdf')
