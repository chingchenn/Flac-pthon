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
from scipy import interpolate

#-----------------------------------SETTING------------------------------------
padding = 30
dx = 0.5
dz = 0.5
# Model size
xmin = 0
xmax = 1200
zmin = -300
zmax = 0
#------------------------------------------------------------------------------
model = sys.argv[1]
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec

for i in range(1,end+1):
    frame = i
    grid_x, grid_z = fd.make_grid(xmin-padding, xmax+padding, zmin-padding, zmax+padding, dx, dz)
    x, z= fl.read_mesh(frame)
    ele_x,ele_z=flac.elem_coord(x, z)
    points = np.vstack((ele_x.flat, ele_z.flat)).T
    
    values = fl.read_visc(frame)
    #vis = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
    vis = fd.gaussian_interpolation2d(ele_x, ele_z, values, grid_x, grid_z)
    f = fd.clip_topo(grid_x, grid_z, vis, x, z)
    fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_visc','/home/jiching/geoflac/data/',
        grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])
    
    mx, mz, mage, mphase, idm, a1, a2, ntriag = fl.read_markers(frame)
    points = np.vstack((mx.flat, mz.flat)).T
    f0 = interpolate.griddata(points, mphase.flatten(), (grid_x, grid_z), method='nearest')
    f0=f0.astype(np.float32)
    f = fd.clip_topo(grid_x, grid_z, f0, x, z)
    fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_ph','/home/jiching/geoflac/data/',
        grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f0[~np.isnan(f)])
    
    # points = np.vstack((x.flat, z.flat)).T
    # values = fl.read_temperature(frame)
    # temp = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
    # f = fd.clip_topo(grid_x, grid_z, temp, x, z)
    # fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_temperature','/home/jiching/geoflac/data/',
    #     grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])
