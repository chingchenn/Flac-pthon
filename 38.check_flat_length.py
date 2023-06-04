#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 4 13:11:24 2023
@author: jiching
"""
import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as fd
import matplotlib.pyplot as plt
from numpy import unravel_index

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/home/jiching/test_geoflac/geoflac/'
#path = '/home/jiching/geoflac_T/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/04model/'
#path = '/scratch2/jiching/'
#path = 'F:/model/'
#savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
#figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
model = sys.argv[1]
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
bwith = 3

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
              "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
              "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
              "#F67280","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
cm = plt.cm.get_cmap('RdYlBu_r')

def trench(end=end):
    trench_x=np.zeros(end)
    trench_z=np.zeros(end)
    trench_index=np.zeros(end)
    arc_x=np.zeros(end)
    arc_z=np.zeros(end)
    arc_index=np.zeros(end)
    for i in range(1,end):
        x,z = fl.read_mesh(i)
        sx,sz=fd.get_topo(x,z)
        arc_ind,trench_ind=fd.find_trench_index(z)
        trench_index[i]=trench_ind
        trench_x[i]=sx[trench_ind]
        trench_z[i]=sz[trench_ind]
        arc_index[i]=arc_ind
        arc_x[i]=sx[arc_ind]
        arc_z[i]=sz[arc_ind]
    return trench_index,trench_x,trench_z,arc_index,arc_x,arc_z
trench_index,trench_x,trench_z,arc_index,arc_x,arc_z = trench()
#def flat_slab_duration():
phase_oceanic = 3;phase_ecolgite = 13
bet = 2;find_flat_dz1=[];find_flat_dz2=[];flat_slab_length=[];flat_slab_depth=[];flat_time=[];flat_length=[]; flat_depth=[]
for i in range(1,end): # 1. find oceanic crust element in each time step 
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
    x_ocean = mx[((phase==phase_ecolgite)+(phase==phase_oceanic))*(mz>-300)]
    z_ocean = mz[((phase==phase_ecolgite)+(phase==phase_oceanic))*(mz>-300)]
    if trench_z[i]> -2 or min(z_ocean)>-200:
        continue
    start = math.floor(trench_x[i]-50)
    final = math.floor(np.max(x_ocean))
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    if len(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])==0:
        continue
    kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
    x_ocean = x_ocean[z_ocean<kk]
    z_ocean = z_ocean[z_ocean<kk]
    for yy,xx in enumerate(x_grid):
        if len(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])==0:
            continue
        oz[yy] = np.average(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        px = xx
    oxx=ox[ox>start]
    oz=oz[ox>start]
    ox=oxx
### =========================== polynomial ===========================


    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(i)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    ax.pcolormesh(ele_x,-ele_z,phase,cmap=phase15,vmin=1, vmax=20)
# # 2. find the fitting of 4th order polynimial of oceanic crust in each time step 
    z1=np.polyfit(ox,oz,5) # 4th order
    p4=np.poly1d(z1)
    w1=p4(ox)
    p3=np.polyder(p4,1) # find f'(x)
    p2=np.polyder(p4,2) # find f"(x)
    w2=p3(ox)
    w3=p2(ox)
    cc=-1;ff1=[];ff1z=[]
    for rr,oo in enumerate(w2): # find slope < 0.2 
        if cc*(oo+0.2)<0:
            ff1.append(ox[rr])
            ff1z.append(oz[rr])
            ax.scatter(ox[rr],-oz[rr],s=50,color='w')
        cc = oo+0.2
    if len(ff1)>=2 and (ff1[-1]-ff1[-2])>10 and ff1[-2]>400:
        flat_time.append(fl.time[i])
        flat_length.append(ff1[-1]-ff1[-2])
        flat_depth.append(np.average(w1[(ox>=ff1[-2])*(ox<ff1[-1])]))
    mm=-1;ff2=[]
    for pp,uu in enumerate(w3): # find inflection points
        if mm*uu<0:
            ff2.append(ox[pp])
        mm = uu  
    if len(ff2)>1 and (ff2[1]-ff2[0])>80 and ff2[0]>start:
        find_flat_dz2.append(fl.time[i])
        if len(ff1)>2 and (ff1[-1]-ff1[-2])>50:
            find_flat_dz1.append(fl.time[i])
            flat_slab_length.append(ff1[-1]-ff1[-2])
            depth=np.average(w1[(ox>=ff1[-2])*(ox<ff1[-1])])
            flat_slab_depth.append(depth)
    


    ax.set_ylim(300,-0)
    ax.set_xlim(200,800)
    ax.scatter(ox,-w1,s=50,color='k')
    ax.scatter(ff1[-1],-ff1z[-1],s=80,color='r')
    ax.scatter(ff1[-2],-ff1z[-2],s=80,color='r')
    ax.set_aspect('equal')
    ax.set_title(fl.time[i],fontsize = 30)
    # xt,zt = fl.read_mesh(frame)
    # temp = fl.read_temperature(frame)
    # ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
    fig.savefig(figpath+'check_slab_length_of_'+str(i)+'.png')
    #return find_flat_dz2,find_flat_dz1,flat_slab_length,flat_slab_depth,flat_time,flat_length,flat_depth
name=model+'_flat_time_len_check'
fs.save_3txt(name,savepath,flat_time,flat_length,flat_depth)