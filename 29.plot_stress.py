#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:48:01 2022

@author: grace
"""


import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
from Main_creat_database import oceanic_slab,nodes_to_elements

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/'
#path = 'F:/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
model = sys.argv[1]
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
bwith = 3

frame = 125
stressxx = fl.read_sxx(frame)
x,z = fl.read_mesh(frame)
ele_x, ele_z = nodes_to_elements(x,z)
fig2, (ax) = plt.subplots(1,1,figsize=(20,22))
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378',
    '#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D',
    '#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
for kk,depth in enumerate([10,20,30,40,50,60,70]):
    sxx = stressxx[:,depth] * 1e2 # kb to MPa
    ax.scatter(ele_x[:,depth],sxx,c = newcolors[kk], lw = 4, label = str(depth*2)+' km')
ax.legend(fontsize = 16)
ax.set_xlim(550,850)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.grid()
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
fig2.savefig(figpath+model+'_sxx_in_depth.png')
