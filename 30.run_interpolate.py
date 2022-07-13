#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:48:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import function_for_flac as fd
import matplotlib.pyplot as plt
from scipy import interpolate


#-------------------------------------------------------------------
model = sys.argv[1]
frame = int(sys.argv[2])
#model = 'cH1404'
#frame = 120
plt.rcParams["font.family"] = "Times New Roman"
path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Volumes/SSD500/data/'
figpath='/home/jiching/geoflac/figure/'
#figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'
os.chdir(path+model)

fl = flac.Flac()
frame = 120
grid_x, grid_z = fd.make_grid(0, 1200, 0, 300, 0.5, 0.5)
values = fl.read_sxx(frame)
x, z= fl.read_mesh(frame)
ele_x,ele_z=flac.elem_coord(x, z)
points = np.vstack((ele_x.flat, ele_z.flat)).T
grid_z1 = interpolate.griddata(points, values, (grid_x, grid_z), method='linear')