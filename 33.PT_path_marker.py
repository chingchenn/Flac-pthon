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
Cocos = 0
Nazca = 1
if Cocos:
    model = 'Cocos_a0906'
    marker_list = [186769,186775,704604]#,928985,960783,973829,992876,928808,928935,889439,882174,886460,864250,861514,859346]
    xmin,xmax=500,950
if Nazca: # 52238,51380,501242, 48622, 60366, 61268, 63076, 62169, 62169, 58564, 59471, 57672, 54078, 52240, 52225, 422137, 732540
    model = 'Ref_Nazca'
    marker_list = [52238,51380,501242, 48622, 60366, 61268, 63076, 62169, 62169, 58564, 59471, 57672, 54078, 52240, 52225, 422137, ]
    xmin,xmax=300,500
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
#figpath='/Users/ji-chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
'''
Cocos9
[248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900]
'''
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
phase15= matplotlib.colors.ListedColormap(phase_color)
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
x = np.linspace(710,1050)
y = -1.25/350*x+5
ax2.plot(x,y,c='#FF9900',lw=5) # solidus

x = np.linspace(680,1050)
y = 0.65/400*x-0.45625
ax2.plot(x,y,c='#FF9900',lw=5) # solidus
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
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)

for kk,idd in enumerate(marker_list):
    px = np.zeros(end)
    pz = np.zeros(end)
    pPre = np.zeros(end)
    pTemp = np.zeros(end)
    pmark= np.zeros(end)
    for ii,qq in enumerate(range(21,end,5)):
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(qq)
        temp = fl.read_temperature(qq-20)
        if len(mx[ID == idd])==0:
            continue
        px[ii] = mx[ID == idd]
        pz[ii] = mz[ID == idd]
        pmark[ii] = phase[ID == idd]
        pPre[ii] = 3300*10*(-pz[ii])*1000/1e9
        pTemp[ii] = flac.marker_interpolate_node(ntriag[ID==idd], a1[ID==idd], a2[ID==idd], fl.nz, temp)
        pTemp[ii] = pTemp[ii] - 0.399*pz[ii]
    # ax2.scatter(pTemp[pTemp>0],pPre[pTemp>0],c=pmark[pTemp>0],cmap = phase19,s = 40,vmin = 1,vmax = 20)
    axdep.scatter(pTemp[pTemp>0],-pz[pTemp>0],c=pmark[pTemp>0],cmap = phase15,s = 40,vmin = 1,vmax = 20)
    # ax2.set_title(str(idd),fontsize = labelsize)
    
    # fig,ax = plt.subplots(1,1,figsize=(10,10))
    # ax.scatter(px[px>0],pz[px>0],color = colors[kk],s = 30)
    # ax.set_xlim(xmin,xmax)
    # ax.set_ylim(-150,0)
    # ax.set_aspect('equal')
    # ax.grid()
    # ax.set_title(str(idd),fontsize = labelsize)
    # ax.tick_params(axis='x', labelsize=labelsize)
    # ax.tick_params(axis='y', labelsize=labelsize)

