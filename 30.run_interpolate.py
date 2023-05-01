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
zmax = 30
#------------------------------------------------------------------------------
model = sys.argv[1]
frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
path = '/scratch2/jiching/04model/'
savepath='/home/jiching/geoflac/data/'
path = '/scratch2/jiching/data/'
figpath='/home/jiching/geoflac/figure/'
path = '/scratch2/jiching/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
#------------------------------------------------------------------------------
run_interpolatation = 1
phasein     = 1
tin         = 0
vis         = 0
gravity     = 0
if not os.path.isdir(path+model+'/phase_vis'):
    os.mkdir(path+model+'/phase_vis')
#-----------------------------------interpolatation----------------------------
if run_interpolatation:
    for i in range(frame,frame+1):
        frame = i
        grid_x, grid_z = fd.make_grid(xmin-padding, xmax+padding, zmin-padding, zmax+padding, dx, dz)
        x, z= fl.read_mesh(frame)
        ele_x,ele_z=flac.elem_coord(x, z)
        points = np.vstack((ele_x.flat, ele_z.flat)).T
        
        values = fl.read_visc(frame)
        vis = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
        #vis = fd.gaussian_interpolation2d(ele_x, ele_z, values, grid_x, grid_z)
        f = fd.clip_topo(grid_x, grid_z, vis, x, z)
        fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_visc_linear','/scratch2/jiching/data/',
            grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])
        
        mx, mz, mage, mphase, idm, a1, a2, ntriag = fl.read_markers(frame)
        points = np.vstack((mx.flat, mz.flat)).T
        f0 = interpolate.griddata(points, mphase.flatten(), (grid_x, grid_z), method='nearest')
        f0 = f0.astype(np.float32)
        f = fd.clip_topo(grid_x, grid_z, f0, x, z)
        fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_ph','/scratch2/jiching/data/',
            grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f0[~np.isnan(f)])
        
        # points = np.vstack((x.flat, z.flat)).T
        # values = fl.read_temperature(frame)
        # temp = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
        # f = fd.clip_topo(grid_x, grid_z, temp, x, z)
        # fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_temperature','/home/jiching/geoflac/data/',
        #     grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])

#---------------------------------------GMT plotting--------------------
if phasein:
    phase_interpolate_file = savepath+model+'_frame_'+str(frame)+'interpolate_ph.txt'
    cmd = '''
cp /home/jiching/GMT/phase19.cpt .
tail -n +2 %(phase_interpolate_file)s | xyz2grd -Gphase_vis/%(model)s_frame%(frame)s_phase.grd -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f

gmt set FONT_ANNOT_PRIMARY          8p,4,#0D057A \
        FONT_TITLE                  20p,4,#0D057A 

gmt begin %(figpath)s/%(model)s_frame%(frame)s_phase pdf
gmt grdimage -JX9c/3c -R0/1200/-300/0 phase_vis/%(model)s_frame%(frame)s_phase.grd -Cphase19.cpt -BWSne -Bxa100f50+l"Distance (Myr)" -Bya100+l"Depth (km)"
#gmt grdcontour temp.grd -Cjet.cpt -A200 -W1p -R0/1200/-300/0
gmt end

rm phase_vis/%(model)s_frame%(frame)s_phase.grd
''' %locals()
    #print(cmd)
    os.system(cmd)
