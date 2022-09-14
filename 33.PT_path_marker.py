#!/usr/bin/env python
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
bwith = 3
colors = ["#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#F67280","#00FF00","#FFFF00","#7158FF"]

phase_color= ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(phase_color)
### 248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900
mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(frame)
px = np.zeros(end)
pz = np.zeros(end)
pres = -fl.read_pres(frame)* 0.1  #kbar to GPa
temp = fl.read_temperature(frame)
den = fl.read_density(frame)
for kk,idd in enumerate([248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900]):
# for kk,idd in enumerate([248644,617652, 485634]):
    px = np.zeros(end)
    pz = np.zeros(end)
    pPre = np.zeros(end)
    pTemp = np.zeros(end)
    pmark= np.zeros(end)
    for ii,qq in enumerate(range(1,120)):
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(qq)
        if len(mx[ID == idd])==0:
            continue
        px[ii] = mx[ID == idd]
        pz[ii] = mz[ID == idd]
        pmark[ii] = phase[ID == idd]
        # pPre[ii] = flac.marker_interpolate_elem(ntriag[ID==idd],fl.nz,pres)
        pPre[ii] = 3000*10*(-pz[ii])*1000/1e9
        pTemp[ii] = flac.marker_interpolate_node(ntriag[ID==idd], a1[ID==idd], a2[ID==idd], fl.nz, temp)
        pTemp[ii] = pTemp[ii] - 0.399*pz[ii]
    
    # fig,ax = plt.subplots(1,1,figsize=(10,10))
    # ax.scatter(px[px>0],pz[px>0],color = colors[kk],s = 30)
    # ax.set_xlim(700,950)
    # ax.set_ylim(-150,0)
    # ax.grid()
    # ax.set_title(str(idd),fontsize = labelsize)
    # ax.tick_params(axis='x', labelsize=labelsize)
    # ax.tick_params(axis='y', labelsize=labelsize)
    fig2,ax2 = plt.subplots(1,1,figsize=(10,10))
    
    x = np.linspace(0,514)
    y = -0.0375 * x + 20.1
    basalt_change = (50/255, 200/255, 180/255)
    ax2.plot(x,y,c=basalt_change,lw=8)
    x = np.linspace(515,1300)
    y = 0.0022 * x - 0.3
    ax2.plot(x,y,c=basalt_change,lw=8,label='basalt-eclogite')
    pressure_limit = 4 # GPa
    pressure=np.linspace(0,pressure_limit,100)
    sss=np.zeros(len(pressure))
    for q,dd in enumerate(pressure):
        if dd<1:
            ss=1050-420*(1-np.exp(-dd*3.3))
        elif dd>2.38:
            ss=(dd+14)*43
        else:
            ss=630+26*((-dd)**2)/2 
        sss[q] = ss
    ax2.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')
    ax2.set_ylim(0,pressure_limit)
    ax2.set_xlim(0,1200)
    axdep = ax2.twinx()
    axdep.set_ylim(0,pressure_limit*1e9/3300/10/1e3)
    ax2.legend(fontsize=labelsize)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax2.tick_params(axis='x', labelsize=labelsize)
    ax2.tick_params(axis='y', labelsize=labelsize)
    ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=labelsize)
    ax2.set_ylabel('Pressure (GPa)',fontsize=labelsize)
    axdep.set_ylabel('Depth (km)',fontsize=labelsize)
    # ax2.text(80,3.5,'Basalt',fontsize=36)
    # ax2.text(570,5.5,'Eclogite',fontsize=36)
    # ax.text(480,2.5,'Phase boundary from',fontsize=fontsize)
    # ax.text(520,2.1,'Hacker et al. (2003)',fontsize=fontsize)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    
    ax2.scatter(pTemp[pTemp>0],pPre[pTemp>0],c=pmark[pTemp>0],cmap = phase19,s = 40,vmin = 1,vmax = 20)
    ax2.grid()
    ax2.set_title(str(idd),fontsize = labelsize)
    ax2.tick_params(axis='x', labelsize=labelsize)
    ax2.tick_params(axis='y', labelsize=labelsize)
    

    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.scatter(px[px>0],pz[px>0],color = colors[kk],s = 30)
    ax.set_xlim(700,950)
    ax.set_ylim(-150,0)
    ax.grid()
    ax.set_title(str(idd),fontsize = labelsize)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    fig.close()
