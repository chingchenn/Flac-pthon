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
model = sys.argv[1]
frame = int(sys.argv[2])
padding = 30
dx = 0.5
dz = 0.6
# Model size
if model=='Ref_Cocos':
    xmin = 500
    xmax = 900
    zmin = -150
    zmax = 20
    path='/home/jiching/geoflac/'
if model=='Nazca_a0702':
    xmin = 250
    xmax = 1000
    zmin = -200
    zmax = 20
    path = '/scratch2/jiching/04model/'
#------------------------------------------------------------------------------
savepath='/home/jiching/geoflac/data/'
savepath = '/scratch2/jiching/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/scratch2/jiching/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
#------------------------------------------------------------------------------
run_interpolatation = 1
phasein     = 0
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
        
        #### melting points ####
        melt = fl.read_fmelt(frame)  
        magma = fl.read_fmagma(frame)
        meltpoint_x = ele_x[melt!=0]
        meltpoint_z = ele_z[melt!=0]
        meltpoint = melt[melt!=0]
        fs.save_3txt(model+'_frame_'+str(frame)+'melting_points',savepath,
                meltpoint_x,meltpoint_z,meltpoint)

        ### END ###
        values = fl.read_visc(frame)
        vis = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
        #vis = fd.gaussian_interpolation2d(ele_x, ele_z, values, grid_x, grid_z)
        f = fd.clip_topo(grid_x, grid_z, vis, x, z)
        fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_visc_linear',savepath,
            grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])
        
        mx, mz, mage, mphase, idm, a1, a2, ntriag = fl.read_markers(frame)
        for kk in range(len(mx)): 
            if mz[kk] > -5 and mphase[kk]==14: 
                mphase[kk]=1 
        points = np.vstack((mx.flat, mz.flat)).T
        f0 = interpolate.griddata(points, mphase.flatten(), (grid_x, grid_z), method='nearest')
        f0 = f0.astype(np.float32)
        f = fd.clip_topo(grid_x, grid_z, f0, x, z)
        fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_ph',savepath,
            grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f0[~np.isnan(f)])
        
        points = np.vstack((x.flat, z.flat)).T
        values = fl.read_temperature(frame)
        temp = interpolate.griddata(points, values.flatten(), (grid_x, grid_z), method='linear')
        #f = fd.clip_topo(grid_x, grid_z, temp, x, z)
       # fs.save_3txt(model+'_frame_'+str(frame)+'interpolate_temperature',savepath,
       #     grid_x[~np.isnan(f)],grid_z[~np.isnan(f)],f[~np.isnan(f)])
        fs.save_3txt(model+'_frame_'+str(frame)+'temperature',savepath,
            x.flatten(),z.flatten(),values.flatten())
        #f = open(model+'_frame_'+str(frame)+'temperature.txt', 'w')
        #f.write('%d %d\n' % x.shape)
        #flac.printing(x, z, values, stream=f)
        #f.close()

#---------------------------------------GMT plotting--------------------
if phasein:
    phase_interpolate_file = savepath+model+'_frame_'+str(frame)+'interpolate_ph.txt'
    temperature_file = savepath+model+'_frame_'+str(frame)+'interpolate_temperature.txt'
    melting_file = savepath+model+'_frame_'+str(frame)+'melting_points.txt'
    time=str(int(frame*0.2))
    cmd = '''
cp /scratch2/jiching/GMT/phase8.cpt .
cp /scratch2/jiching/GMT/cw.cpt .
tail -n +2 %(phase_interpolate_file)s | xyz2grd -Gphase_vis/%(model)s_frame%(frame)s_phase.grd -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f
tail -n +2 %(temperature_file)s | surface -Gphase_vis/%(model)s_frame%(frame)s_temp.grd -Ll0 -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f

gmt set FONT_ANNOT_PRIMARY          8p,4,#000000 \
        FONT_TITLE                  14p,4,#000000 

gmt begin %(figpath)s%(model)s_frame%(frame)s_phase pdf
gmt grdimage -JX15c/5c -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f phase_vis/%(model)s_frame%(frame)s_phase.grd -Cphase8.cpt -BWSne+t"Time %(time)s Myr" -Bxa100f50+l"Distance (km)" -Bya50+l"Depth (km)"
#gmt basemap -JX18c/6c -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -BWSne+t"Time %(time)s Myr" -Bxa100f50+l"Distance (km)" -Bya50+l"Depth (km)"
#cat %(phase_interpolate_file)s | gmt plot -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -Cphase8.cpt -Sc0.05
gmt grdcontour phase_vis/%(model)s_frame%(frame)s_temp.grd -Ccw.cpt -A200 -Wa1p -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f
gmt end

gmt begin %(figpath)s%(model)s_frame%(frame)s_phase png
gmt grdimage -JX18c/6c -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f phase_vis/%(model)s_frame%(frame)s_phase.grd -Cphase8.cpt -BWSne+t"Time %(time)s Myr" -Bxa100f50+l"Distance (km)" -Bya50+l"Depth (km)"
cat %(melting_file)s | gmt plot -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f -G220/20/60 -Sc0.1
gmt end
pwd
mv phase_vis/%(model)s_frame%(frame)s_phase.grd %(savepath)s%(model)s_frame%(frame)s_phase.grd
''' %locals()
    print(cmd)
    os.system(cmd)
