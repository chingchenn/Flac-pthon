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
import function_savedata as fs
import matplotlib.pyplot as plt
from scipy import interpolate


#-------------------------------------------------------------------
model = sys.argv[1]
#frame = int(sys.argv[2])
frame = 120
path='/home/jiching/geoflac/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)

fl = flac.Flac()


padding = 30
grid_x, grid_z = fd.make_grid(0-padding, 1200+padding, -300-padding, 0+padding, 0.5, 0.5)
values = fl.read_visc(frame)
x, z= fl.read_mesh(frame)
ele_x,ele_z=flac.elem_coord(x, z)
points = np.vstack((ele_x.flat, ele_z.flat)).T
grid_z1 = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_visc','/home/jiching/geoflac/data/',
    grid_x[~np.isnan(grid_z1)],grid_z[~np.isnan(grid_z1)],grid_z1[~np.isnan(grid_z1)])

mx, mz, mage, mphase, idm, a1, a2, ntriag = fl.read_markers(frame)
points = np.vstack((mx.flat, mz.flat)).T
f0 = interpolate.griddata(points, mphase.flatten(), (grid_x, grid_z), method='nearest')
f0=f0.astype(np.float32)
f = fd.clip_topo(grid_x, grid_z, f0, x, z)
fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_ph','/home/jiching/geoflac/data/',
    grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f0[~np.isnan(f)])

