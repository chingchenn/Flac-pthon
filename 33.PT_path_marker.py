# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 10:13:28 2022

@author: grace
"""

import flac
import math
import sys, os
import matplotlib
import numpy as np
from matplotlib import cm
from scipy import interpolate
import function_for_flac as fd
# import Main_creat_database as Md
import function_savedata as fs
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
model = 'Cocos9'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
#figpath='/Users/ji-chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1

frame = 6
labelsize = 20
colors = ["#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#F67280","#00FF00","#FFFF00","#7158FF"]
### 248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900
mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(frame)
px = np.zeros(end)
pz = np.zeros(end)
for kk,idd in enumerate([248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900]):
    px = np.zeros(end)
    pz = np.zeros(end)
    for ii,qq in enumerate(range(1,120)):
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(qq)
        if len(mx[ID == idd])==0:
            continue
        px[ii] = mx[ID == idd]
        pz[ii] = mz[ID == idd]
    
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.scatter(px[px>0],pz[px>0],color = colors[kk],s = 30)
    ax.set_xlim(700,950)
    ax.set_ylim(-150,0)
    ax.grid()
    ax.set_title(str(idd),fontsize = labelsize)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    fig.close()
    

