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
# plt.figure(figsize=(9, 1.5))
fig, ax = plt.subplots(figsize=(2, 17))
# fig.subplots_adjust(bottom=0.5)

cmap = mpl.cm.rainbow
# norm = mpl.colors.Normalize(vmin=0, vmax=200)

# cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
#                                 norm=norm,spacing='uniform',
#                                 orientation='vertical')
# # cb1.set_title('Pa s',fontsize=20)
# cb1.set_label('Some Units')
# fig.show()

# fig, ax = plt.subplots(figsize=(1, 20))
# fig.subplots_adjust(bottom=0.5)

# cmap = mpl.cm.jet
bounds =[0,200,400,600,800,1000,1200]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb2 = mpl.colorbar.ColorbarBase(ax,cmap=cmap,
                                 norm=norm,
                                orientation='vertical')
cb2.ax.tick_params(labelsize=40)
# # cb2.set_label("Discrete intervals with extend='both' keyword")
# fig.show()

# # plt.subplots(figsize=(1, 3))
# ax0 = plt.axes([0.3, 0.1, 0.1, 0.8])

# # add colorbars
# cmap = mpl.cm.jet
# # normalize = mpl.colors.Normalize(0, 10)
# boundarynorm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
# cb0 = mpl.colorbar.ColorbarBase(ax=ax0, norm=boundarynorm,cmap=cmap)


# # add labels
# cb0.set_label('BoundaryNorm')
