#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:28:16 2021

@author: jiching
"""
import math
import flac
import os,sys
import numpy as np
from scipy import interpolate
# import matplotlib
# matplotlib.use('Agg')
import function_savedata as fs
import matplotlib.pyplot as plt
import function_for_flac as f2
#=========================setting=============================
model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model+'/'
#path = 'D:/model/'+model+'/'
#path = '/scratch2/jiching/sem02model/'+model+'/'
path = '/scratch2/jiching/03model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
#=========================Parameters=========================
phase_oceanic = 3
phase_ecolgite = 13
bet = 2
find_flat_dz1=[]
find_flat_dz2=[]
flat_slab_length=[]
flat_slab_depth=[]
figg2=0
#=========================main code===========================
for i in range(1,end):
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
    ## In this code, we considered the marker phase, not the element phase
    trench_ind = np.argmin(z[:,0]) 
    x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]
    x_ocean = mx[(phase==phase_ecolgite)+(phase==phase_oceanic)]
    z_ocean = mz[(phase==phase_ecolgite)+(phase==phase_oceanic)]
    if z_trench> -2 or min(z_ocean)>-200:
        continue
    start = math.floor(x_trench)
    final = math.floor(np.max(x_ocean))
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    #find initial basalt depth to remove the weage basalt
    kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
    x_ocean = x_ocean[z_ocean<kk]
    z_ocean = z_ocean[z_ocean<kk]
    # interplate to the grid length "bet"
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
    z1=np.polyfit(ox,oz,4)
    p4=np.poly1d(z1)
    w1=p4(ox)
    p3=np.polyder(p4,1)
    p2=np.polyder(p4,2)
    w2=p3(ox)
    w3=p2(ox)
    if figg2:
        ### This figure show the origin point fitting with polynomail function
        fig2,(q1,q2,q3)= plt.subplots(3,1,figsize=(9,12))
        q1.plot(ox,w1,c='k',lw=3,label ='4st polynomial')
        q1.scatter(ox,oz,c='cyan',s=20,label ='average marker points')           
        q2.plot(ox,w2,c='k',label ='3st polynomial')
        q2.axhline(y=-0.1, xmin=0, xmax=800,color='r',linestyle='dashed')
        q3.plot(ox,w3,c='k',label ='2st polynomial')
        q3.plot([start,final],[0,0],'--',zorder=0,color='red')
        q1.set_title('frame='+str(i))
        q2.set_xlim(start,final)
        q1.set_xlim(start,final) 
        q3.set_xlim(start,final) 
        q1.grid();q2.grid();q3.grid() 
        q1.tick_params(axis='x', labelsize=16)
        q1.tick_params(axis='y', labelsize=16)
        q2.tick_params(axis='x', labelsize=16)
        q2.tick_params(axis='y', labelsize=16)
        q3.tick_params(axis='y', labelsize=16)
        q3.set_xlabel('Distance (km)')
        fig2.savefig('/home/jiching/geoflac/figure/'+model+'_frame='+str(i)+'_fig2.png')
    cc=-1;ff1=[]
    for rr,oo in enumerate(w2):
        if cc*(oo+0.1)<0:
            ff1.append(ox[rr])
        cc = oo+0.1
    mm=-1;ff2=[]
    for pp,uu in enumerate(w3):
        if mm*uu<0:
            ff2.append(ox[pp])
        mm = uu  
    if len(ff2)>1 and (ff2[1]-ff2[0])>100 and ff2[0]>start:
        find_flat_dz2.append(fl.time[i])
        if len(ff1)>3 and (ff1[-1]-ff1[-2])>50:
            find_flat_dz1.append(fl.time[i])
            flat_slab_length.append(ff1[-1]-ff1[-2])
            depth=np.average(w1[(ox>=ff1[-2])*(ox<ff1[-1])])
            flat_slab_depth.append(depth)
fs.save_1txt(str(model)+'_flatslab_duration2','/home/jiching/geoflac/data',find_flat_dz2)
fs.save_3txt(str(model)+'_flatslab_time_len','/home/jiching/geoflac/data',find_flat_dz1,flat_slab_length,flat_slab_depth)
