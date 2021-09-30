#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:14 2021

@author: ji-chingchen
"""


import flac
import os
import numpy as np
import matplotlib.pyplot as plt

# model = 'w1202'
# #path = '/home/jiching/geoflac/'+model+'/'
# path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
# path = '/Volumes/My Book/model/'+model+'/'
# os.chdir(path)
# # fl = flac.Flac()
# fl = flac.Flac();end = fl.nrec
# total_xx = np.zeros(end)
# for kk in range(1,end):
#     xx = fl.read_sxx(kk)
#     x,z = fl.read_mesh(kk)
#     zz = z[0,:]
#     bxx = xx[0,:]
#     for yy in range(len(bxx)):
#         # print(yy,zz[yy],zz[yy+1])
#         dl = abs(zz[yy]-zz[yy+1])
#         total_xx[kk] += bxx[yy] * dl
#         # print(total_xx[kk], bxx[yy])
# fig, ax= plt.subplots(1,1,figsize=(12,8))        
# ax.plot(fl.time,total_xx,c='k')
# ax.set_xlim(0,fl.time[-1])
# ax.tick_params(axis='x', labelsize=16)
# ax.tick_params(axis='y', labelsize=16)
model='s0209'
path = '/Users/ji-chingchen/Desktop/model/data/forc'+model
temp1=np.loadtxt(path+".txt")
time,forc,forc_r = temp1.T
fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
fig, (ax)= plt.subplots(1,1,figsize=(12,8))     
ax.scatter(time,forc)
x_grid = np.arange(0.2,time[-1],0.1)
ox = np.zeros(len(x_grid))
oz = np.zeros(len(x_grid))
px = 0
for yy,xx in enumerate(x_grid):
    oz[yy] = np.average(forc[(time>=px)*(time<xx)])
    ox[yy] = np.average(time[(time>=px)*(time<xx)])
    # ox[yy] = np.average(x_ocean[(x_ocean>=px) *(x_ocean<=xx)])
    px = xx
ax.plot(ox,oz,c='k')
ax.set_xlim(0,time[-1])
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.plot([time[0],time[-1]],[-6.31*10**12,-6.31*10**12],'r--')
oo = np.zeros(len(ox))
for qq in range(1,len(ox)):
    oo[qq] = (oz[qq]-oz[qq-1]) /(ox[qq]-ox[qq-1])
# ax2.plot(ox,oo,c='k')
# ax2.set_xlim(0,time[-1])
# ax2.set_ylim(-0.5*10**13,0.5*10**13)
# ax2.tick_params(axis='x', labelsize=16)
# ax2.tick_params(axis='y', labelsize=16)
# ax2.plot([time[0],time[-1]],[0,0],'r--')
