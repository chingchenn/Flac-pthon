#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 22:24:33 2023

@author: chingchen
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
from numpy import unravel_index
import matplotlib.colors as mcolors

path = '/home/jiching/geoflac/'
path='/Users/chingchen/Desktop/model/'
#savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
savepath='/Users/chingchen/Desktop/data/'
#figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
figpath='/Users/chingchen/Desktop/figure/'
# model = sys.argv[1]


plt.rcParams["font.family"] = "Helvetica"
bwith = 3

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

model='Ref_Cocos'
#model='Nazca_aa06'
os.chdir(path+model)
fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
# end = 250
time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
rainbow = cm.get_cmap('gray_r',end)
meltcolor = cm.get_cmap('turbo',end)
newcolors = rainbow(np.linspace(0, 1, end))
time_color = meltcolor(np.linspace(0,1,end))


size = 120

fig, (ax1,ax2,ax5)= plt.subplots(3,1,figsize=(15,15),gridspec_kw={'height_ratios':[3,3,3]})
xxx_trench = np.max(x_trench)


#### time 
i = 100     ## 10 Myr 
i2 = 100   ## 20 Myr 
i3 = 200   ## 40 Myr



def x2dis1(x):
    return x -x_trench[i]
def dis2x1(x):
    return x +x_trench[i]
# ax3 = ax1.secondary_xaxis('top', functions=(dis2x1,x2dis1))
def x2dis2(x):
    return x -x_trench[i2]
def x2dis3(x):
    return x -x_trench[i3]



### eclogite    
x, z = fl.read_mesh(i)
ele_x, ele_z = flac.elem_coord(x,z)
phase = fl.read_phase(i)

for xx in range(len(ele_x)):
    if True in (phase[xx,:]==phase_eclogite):
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_eclogite):
                x_change_ec=ele_x[xx,zz]
                z_change_ec=-ele_z[xx,zz]
                break
        break
ax1.scatter(x2dis1(x_change_ec),z_change_ec,c='darkblue',s=size,label = 'basalt')


### sediment
for xx in range(len(ele_x)):
    print(xx)
    if True in (phase[xx,:]==phase_schist):
        print(xx)
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_schist):
                x_change_sd=ele_x[xx,zz]
                z_change_sd=-ele_z[xx,zz]
                break
        break
ax1.scatter(x2dis1(x_change_sd),z_change_sd,c='orange',s=size,label = 'sediment')


### serpentinite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_serpentinite):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_serpentinite):
                x_change_serp=ele_x[xx,zz]
                z_change_serp=-ele_z[xx,zz]
                break
        break
ax1.scatter(x2dis1(x_change_serp),z_change_serp,c='purple',s=size,label = 'serpentinite')

### chlorite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_hydratedmantle):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_hydratedmantle):
                x_change_cho=ele_x[xx,zz]
                z_change_cho=-ele_z[xx,zz]
                break
        break
ax1.scatter(x2dis1(x_change_cho),z_change_cho,c='darkgreen',s=size,label = 'chlorite')

####------------------ time 2 ------------------------
### eclogite    
x, z = fl.read_mesh(i2)
ele_x, ele_z = flac.elem_coord(x,z)
phase = fl.read_phase(i2)

for xx in range(len(ele_x)):
    if True in (phase[xx,:]==phase_eclogite):
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_eclogite):
                x_change_ec=ele_x[xx,zz]
                z_change_ec=-ele_z[xx,zz]
                break
        break
ax2.scatter(x2dis2(x_change_ec),z_change_ec,c='darkblue',s=size)


### sediment
for xx in range(len(ele_x)):
    if True in (phase[xx,:]==phase_schist):
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_schist):
                x_change_sd=ele_x[xx,zz]
                z_change_sd=-ele_z[xx,zz]
                break
        break
ax2.scatter(x2dis2(x_change_sd),z_change_sd,c='orange',s=size)


### serpentinite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_serpentinite):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_serpentinite):
                x_change_serp=ele_x[xx,zz]
                z_change_serp=-ele_z[xx,zz]
                break
        break
ax2.scatter(x2dis2(x_change_serp),z_change_serp,c='purple',s=size)

### chlorite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_hydratedmantle):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_hydratedmantle):
                x_change_cho=ele_x[xx,zz]
                z_change_cho=-ele_z[xx,zz]
                break
        break
ax2.scatter(x2dis2(x_change_cho),z_change_cho,c='darkgreen',s=size)




####------------------ time 3 ------------------------
### eclogite    
x, z = fl.read_mesh(i3)
ele_x, ele_z = flac.elem_coord(x,z)
phase = fl.read_phase(i3)

for xx in range(len(ele_x)):
    if True in (phase[xx,:]==phase_eclogite):
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_eclogite):
                x_change_ec=ele_x[xx,zz]
                z_change_ec=-ele_z[xx,zz]
                break
        break
ax5.scatter(x2dis3(x_change_ec),z_change_ec,c='darkblue',s=size)


### sediment
for xx in range(len(ele_x)):
    if True in (phase[xx,:]==phase_schist):
        for zz in range(len(ele_x[0])):
            if (phase[xx,zz]==phase_schist):
                x_change_sd=ele_x[xx,zz]
                z_change_sd=-ele_z[xx,zz]
                break
        break
ax5.scatter(x2dis3(x_change_sd),z_change_sd,c='orange',s=size)


### serpentinite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_serpentinite):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_serpentinite): 
                x_change_serp=ele_x[xx,zz]
                z_change_serp=-ele_z[xx,zz]
                break
        break
ax5.scatter(x2dis3(x_change_serp),z_change_serp,c='purple',s=size)

### chlorite
for xx in range(len(ele_x)-1,0,-1):
    if True in (phase[xx,:]==phase_hydratedmantle):
        for zz in range(len(ele_x[0])-1,0,-1):
            if (phase[xx,zz]==phase_hydratedmantle):
                x_change_cho=ele_x[xx,zz]
                z_change_cho=-ele_z[xx,zz]
                break
        break
ax5.scatter(x2dis3(x_change_cho),z_change_cho,c='darkgreen',s=size)

for ax in [ax1,ax2,ax5]:
    ax.grid()
    ax.set_xlim(0,600)
    ax.set_ylim(150,0)
    ax.set_aspect('equal')
    ax.set_ylabel('Depth (km)',fontsize=26)
    ax.tick_params(labelsize=26,width = 2)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(bwith)
    ax.set_aspect('equal')
# for ax in [ax3]:#,ax4,ax6]:
#     ax.tick_params(labelsize=26,length = 12,direction='in',width = 2)
#     for axis in ['top','bottom','left','right']:
#         ax.spines[axis].set_linewidth(bwith)
    
frame1 = i*0.2
xx,zz = np.loadtxt(savepath+model+'_'+str(frame1)+'_final_slab.txt').T
xx=xx[zz<0]
zz=zz[zz<0]
zz = fd.moving_window_smooth(zz,3)
ax1.plot(xx+10,-zz,color='k',lw=5)

frame2 = i2*0.2
xx,zz = np.loadtxt(savepath+model+'_'+str(frame2)+'_final_slab.txt').T
xx=xx[zz<0]
zz=zz[zz<0]
zz = fd.moving_window_smooth(zz,3)
ax2.plot(xx+20,-zz,color='k',lw=5)

frame3 = i3*0.2
xx,zz = np.loadtxt(savepath+model+'_'+str(frame3)+'_final_slab.txt').T
xx=xx[zz<0]
zz=zz[zz<0]
zz = fd.moving_window_smooth(zz,3)
ax5.plot(xx,-zz,color='k',lw=5)


ax5.set_title('Location of phase change at '+str(frame3)+' Myr',fontsize=28)
ax5.set_xlabel('Distance from trench (km)',fontsize=28)
ax1.set_title('Location of phase change at '+str(frame1)+' Myr',fontsize=28)
ax2.set_title('Location of phase change at '+str(frame2)+' Myr',fontsize=28)
ax1.legend(fontsize = 25)
#fig.tight_layout()
# fig.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Ref_Nazca/'+'Nazca_a0702_2Dtime_series_7.pdf')

### eclogite    
# x_change_ec = np.zeros(end)
# z_change_ec = np.zeros(end)
# time = np.ones(end)
# for i in range(2,end):
#     x, z = fl.read_mesh(i)
#     ele_x, ele_z = flac.elem_coord(x,z)
#     phase = fl.read_phase(i)
#     time[i] = i * 0.2
#     for xx in range(len(ele_x)):
#         if True in (phase[xx,:]==phase_eclogite):
#             for zz in range(len(ele_x[0])):
#                 if (phase[xx,zz]==phase_eclogite):
#                     x_change_ec[i]=ele_x[xx,zz]
#                     z_change_ec[i]=-ele_z[xx,zz]
#                     break
#             break
# cbtime=ax5.scatter(x_change_ec[x_change_ec>0],z_change_ec[x_change_ec>0],c='darkblue')


# ### sediment
# x_change_sd = np.zeros(end)
# z_change_sd = np.zeros(end)
# time = np.ones(end)
# for i in range(100,end):
#     x, z = fl.read_mesh(i)
#     ele_x, ele_z = flac.elem_coord(x,z)
#     phase = fl.read_phase(i)
#     time[i] = i * 0.2
#     for xx in range(len(ele_x)):
#         if True in (phase[xx,:]==phase_schist):
#             for zz in range(len(ele_x[0])):
#                 if (phase[xx,zz]==phase_schist):
#                     x_change_sd[i]=ele_x[xx,zz]
#                     z_change_sd[i]=-ele_z[xx,zz]
#                     break
#             break
# cbtime=ax5.scatter(x_change_sd[x_change_sd>0],z_change_sd[x_change_sd>0],c='orange')


# ### serpentinite
# x_change_serp = np.zeros(end)
# z_change_serp = np.zeros(end)
# time = np.ones(end)
# for i in range(150,end,20):
#     x, z = fl.read_mesh(i)
#     ele_x, ele_z = flac.elem_coord(x,z)
#     phase = fl.read_phase(i)
#     time[i] = i * 0.2
#     for xx in range(len(ele_x)-1,0,-1):
#         if True in (phase[xx,:]==phase_serpentinite):
#             print('-----  '+str(i)+'  -----')
#             for zz in range(len(ele_x[0])-1,0,-1):
#                 if (phase[xx,zz]==phase_serpentinite):
#                     print(xx,zz)
#                     x_change_serp[i]=ele_x[xx,zz]
#                     z_change_serp[i]=-ele_z[xx,zz]
#                     break
#             break
# cbtime=ax5.scatter(x_change_serp[x_change_serp>0],z_change_serp[x_change_serp>0],c=time[x_change_serp>0],cmap = 'turbo',vmin=0,vmax=40)

# ### chlorite
# x_change_cho = np.zeros(end)
# z_change_cho = np.zeros(end)
# time = np.ones(end)
# for i in range(150,end,20):
#     x, z = fl.read_mesh(i)
#     ele_x, ele_z = flac.elem_coord(x,z)
#     phase = fl.read_phase(i)
#     time[i] = i * 0.2
#     for xx in range(len(ele_x)-1,0,-1):
#         if True in (phase[xx,:]==phase_hydratedmantle):
#             for zz in range(len(ele_x[0])-1,0,-1):
#                 if (phase[xx,zz]==phase_hydratedmantle):
#                     print(xx,zz)
#                     x_change_cho[i]=ele_x[xx,zz]
#                     z_change_cho[i]=-ele_z[xx,zz]
#                     break
#             break
# cbtime=ax5.scatter(x_change_cho[x_change_cho>0],z_change_cho[x_change_cho>0],c=time[x_change_cho>0],cmap = 'turbo',vmin=0,vmax=40)






