#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 15:19:30 2022

@author: chingchen
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

Cocos = 1
Nazca = 0
if Cocos:
    max_dis = 620
    xmin = 500
    xmax = 1100
    model = 'Cocos_x0102'
    marker_list = [186769,186775,704604,928985,960783,973829,992876,928808,928935,889439,882174,886460,864250,861514,859346]

if Nazca:
    max_dis = 515
    xmin = 300
    xmax = 1000
    # model = 'Ref_Nazca'
    model = 'Nazca_a0702'
    marker_list = [52238]#, 51380, 501242, 48622, 60366, 61268, 63076, 62169, 62169, 58564, 59471, 57672, 54078, 52240, 52225, 422137, ]


path='/home/jiching/geoflac/'
path = '/Users/chingchen/Desktop/model/'
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1

phase_uppercrust = 2
phase_oceanic = 3
phase_mantle1 = 4
phase_schist = 5
phase_mantle2 = 8
phase_serpentinite = 9
phase_sediment = 10
phase_sediment_1 = 11
phase_eclogite = 13
phase_lowercrust = 14
phase_hydratedmantle = 16
phase_oceanic_1 = 17
phase_eclogite_1 = 18

'''
Cocos9
[248644, 248647, 247521, 634356, 617652, 485634, 598651, 235136, 573900]
'''

frame = 120
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
def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt


           
fig2,ax2 = plt.subplots(1,1,figsize=(10,10))
x = np.linspace(0,514)
y = -0.0375 * x + 20.1
basalt_change = (50/255, 200/255, 180/255)
ax2.plot(x,y,c=basalt_change,lw=5)
x = np.linspace(515,1300)
y = 0.0022 * x - 0.3
ax2.plot(x,y,c=basalt_change,lw=5,label='basalt-eclogite')
pressure_limit = 5 # GPa
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
ax2.tick_params(labelsize=labelsize)
ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=labelsize)
ax2.set_ylabel('Pressure (GPa)',fontsize=labelsize)
axdep.set_ylabel('Depth (km)',fontsize=labelsize)
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(bwith)

for frame in range(end-2,end):
    # ------ read data from model -----
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    temp = fl.read_temperature(frame)
    temp_ele = temp_elements(temp)
    phase = fl.read_phase(frame)
    pressure = fl.read_pres(frame) * -0.1
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    ind_trench = int(trench_index[frame-1])
    
    # ----- empty array and data -----
    slab_P = np.zeros(len(ele_z))
    slab_T = np.zeros(len(ele_z))
    slab_z = np.zeros(len(ele_z))
    # ----- Start Calculation ------
    for ii,x_ind in enumerate(range(ind_trench-10,len(ele_z))):
        ind_slab = (ele_z[x_ind,:]> -150)*(ele_z[x_ind,:]< -5)*((phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic))
        if not True in ind_slab:
            continue
        if ele_x[x_ind,0]<max_dis:
            top_slab_index = np.where(ind_slab)[0][-1]
        elif ele_x[x_ind,0]<max_dis+410:
            top_slab_index = np.where(ind_slab)[0][-1]
        else:
            continue
        slab_T[x_ind] = temp_ele[x_ind,top_slab_index]
        slab_z[x_ind] = -ele_z[x_ind,top_slab_index]
        slab_P[x_ind] = pressure[x_ind,top_slab_index] # GPa
    slab_z = fd.moving_window_smooth(slab_z, 1)
    axdep.plot(slab_T[slab_T>0],slab_z[slab_T>0],c='green')


# ax2.scatter(slab_T,slab_P,c='#455E45')
for kk,idd in enumerate(marker_list):
    px = np.zeros(end)
    pz = np.zeros(end)
    pPre = np.zeros(end)
    pTemp = np.zeros(end)
    pmark= np.zeros(end)
    for ii,qq in enumerate(range(21,150)):
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(qq)
        temp = fl.read_temperature(qq)
        if len(mx[ID == idd])==0:
            continue
        px[ii] = mx[ID == idd]
        pz[ii] = mz[ID == idd]
        pmark[ii] = phase[ID == idd]
        pPre[ii] = 3300*10*(-pz[ii])*1000/1e9
        pTemp[ii] = flac.marker_interpolate_node(ntriag[ID==idd], a1[ID==idd], a2[ID==idd], fl.nz, temp)
        pTemp[ii] = pTemp[ii] - 0.4*pz[ii]
    # ax2.scatter(pTemp[pTemp>0],pPre[pTemp>0],c=pmark[pTemp>0],cmap = phase15,s = 40,vmin = 1,vmax = 20)
    # axdep.scatter(pTemp[pTemp>0],-pz[pTemp>0],c=pmark[pTemp>0],cmap = phase15,s = 40,vmin = 1,vmax = 20)
    ax2.plot(pTemp[pTemp>0]+0,pPre[pTemp>0],c='purple')
