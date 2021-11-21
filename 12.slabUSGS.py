#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:14:24 2021
@author: jiching
"""
import flac
import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fig = 1
path = '/home/jiching/geoflac/figure/'
#-------------------------------call gmt to cut the trace----------------------------------
cmd = 'cp /scratch2/jiching/GMT/slab/cam_slab2_dep_02.24.18.grd .'
os.system(cmd)
cmd = 'gmt grdtrack -E-100.7/15.5/-97.3/22.0+i0.5k -Gcam_slab2_dep_02.24.18.grd >table.txt'
os.system(cmd)
#------------------------------------------------------------------------------------------
temp2 = np.loadtxt('table.txt')
x,y,z = temp2.T
z1=np.polyfit(x,z,4)
p4=np.poly1d(z1)
w1=p4(x)
#p3=np.polyder(p4,1)
if fig :
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    ax.scatter(x,z,s=2,color='#4169E1')
    ax.plot(x,w1,ls=2,color='k')
    fig.savefig(path+'slab.png')
cmd = 'rm cam_slab2_dep_02.24.18.grd table.txt'
os.system(cmd)
