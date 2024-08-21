#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:15:55 2022

@author: chingchen
"""

#-------------------------------------------------------------------
import sys, os
import numpy as np

import flac
import flac_interpolate as fi
import flac_gravity3 as fg
import matplotlib
import matplotlib.pyplot as plt

###############################################

#model = sys.argv[1]
# frame = int(sys.argv[2])
model = 'Nazca_aa06'
model = 'Nazca_v2_01'
#model = 'Ref_Cocos'
# frame = 10

def make_phase_interploate(model,frame):
    path='/home/jiching/geoflac/'
    path = '/Users/chingchen/Desktop/model/'


    os.chdir(path+model)

    fl = flac.Flac()
    savepath='/scratch2/jiching/data/'
    savepath = '/Users/chingchen/Desktop/data/'
    figpath='/scratch2/jiching/figure/'
    figpath = '/Users/chingchen/Desktop/figure/'
    def find_trench_index(z):
        '''Returns the i index of trench location.'''
        zz = z[:,0]
        # the highest point defines the forearc
        imax = zz.argmax()
        # the trench is the lowest point west of forearc
        i = zz[:imax].argmin()
        return i


    def interpolate_phase(frame, xtrench):
        # domain bounds in km
        fi.xmin = xtrench + left
        fi.xmax = xtrench + right
        fi.zmin = down
        fi.zmax = up

        # resolution in km
        fi.dx = dx
        fi.dz = dz

        xx, zz, ph = fi.interpolate(frame, 'phase')
        return xx, zz, ph
    # domain bounds
    left = -200
    right = 1000
    up = 10
    down = -200
    dx = 1.5
    dz = 0.6
    x, z = fl.read_mesh(frame)
    
    itrench = find_trench_index(z)
    xtrench = x[itrench,0]
    
    # get interpolated phase either from previous run or from original data
    phasefile = 'intp3-phase.%d' % frame
    xx, zz, ph = interpolate_phase(frame, xtrench)
    f = open(phasefile, 'w')
    f.write('%d %d\n' % xx.shape)
    flac.printing(xx, zz, ph, stream=f)
    f.close()
    
    # get topography and gravity at uniform spacing
    px, topo, topomod, gravity = fg.compute_gravity2(frame)
    # convert to km and mGal before saving
    px *= 1e-3
    topo *= 1e-3
    topomod *= 1e-3
    gravity *= 1e5
    gfile = 'topo-grav.%d' % frame
    f = open(gfile, 'w')
    flac.printing(px, topo, gravity, topomod, stream=f)
    f.close()
    
    
    ###############
    model = os.path.split(os.getcwd())[-1]
    phgrd = model+'_phase_%d.grd' % frame
    
    xmin = xtrench + left
    xmax = xtrench + right
    zmin = down
    zmax = up
    aspect_ratio = float(up - down) / (right - left)
    width = 6.5
    height = width * aspect_ratio
    shiftz = height + 0.3
    
    # height of gravity plot
    height2 = 1.0
    
    # gravity grid
    gravgridsize = 50
    gmin = int(gravity.min() / gravgridsize - 1) * gravgridsize
    gmax = int(gravity.max() / gravgridsize + 1) * gravgridsize
    gravann = max(abs(gmin), abs(gmax))
    # topography grid
    topogridsize = 2
    tpmin = int(topo.min() / topogridsize - 1) * topogridsize
    tpmax = int(topo.max() / topogridsize + 1) * topogridsize
    topoann = max(abs(tpmin), abs(tpmax))
    # interval of temperature contours
    cint = 200
    
    if not os.path.exists(phgrd):
        cmd = '''tail -n +2 %(phasefile)s | gmt xyz2grd -G%(phgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f
        mv  %(phgrd)s %(savepath)s
        ''' % locals()
        print(cmd)
        os.system(cmd)


frame_list=[30,60,120,140]
# end = fl.nrec
for frame in range(2,221,2):
#for frame in frame_list:
    make_phase_interploate(model,frame)
    print('---------------end of '+str(frame)+'---------------')
