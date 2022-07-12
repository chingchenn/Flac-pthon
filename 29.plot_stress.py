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
fig2, (ax) = plt.subplots(5,1,figsize=(20,32))
depth1 = 20
depth2 = 80
ax.plot(x[:,depth1],stressxx[:,depth1],c = 'red',lw = 4)
ax.plot(x[:,depth2],stressxx[:,depth2],c = 'b', lw = 4)
fig2.savefig(figpath+model+'_sxx_in_depth.png')
