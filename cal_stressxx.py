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

model = 'w1206'
#path = '/home/jiching/geoflac/'+model+'/'
path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
path = '/Volumes/My Book/model/'+model+'/'
os.chdir(path)
# fl = flac.Flac()
fl = flac.Flac();end = fl.nrec
total_xx = np.zeros(end)
for kk in range(1,end):
    xx = fl.read_sxx(kk)
    x,z = fl.read_mesh(kk)
    zz = z[0,:]
    bxx = xx[0,:]
    for yy in range(len(bxx)):
        # print(yy,zz[yy],zz[yy+1])
        dl = abs(zz[yy]-zz[yy+1])
        total_xx[kk] += bxx[yy] * -dl
        # print(total_xx[kk], bxx[yy])
fig, ax= plt.subplots(1,1,figsize=(12,8))        
ax.plot(fl.time,total_xx,c='k')
ax.set_xlim(0,fl.time[-1])
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)