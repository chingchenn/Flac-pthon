#!/usr/bin/env python

import sys, os
import numpy as np

import flac
import flac_interpolate as fi
import flac_gravity3 as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

#-------------------------------------------------------------------
model = sys.argv[1]
frame = int(sys.argv[2])
plt.rcParams["font.family"] = "Times New Roman"
path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Volumes/SSD500/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'

phasein     = 1
tin         = 0
vis         = 0
gravity     = 0
#-------------------------------------------------------------------
if phasein:
    phase_interpolate_file = savepath+model+'_frame_'+str(frame)+'interpolate_ph.txt'
    cmd = 'tail -n +2 %(phase_interpolate_file)s | xyz2grd -G%(phgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
    print(cmd)
    #os.system(cmd)
