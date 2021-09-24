#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:28:16 2021

@author: jiching
"""
import math
import flac
import os
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as f2
# model = str(sys.argv[1])
# path = '/home/jiching/geoflac/'+model+'/'
model='w1261'
    # path = '/scratch2/jiching/'+model+'/'
    # path = '/home/jiching/geoflac/'+model+'/'
# path = '/Volumes/My Book/model/'+model+'/'
path = '/Volumes/SSD500/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
# flk = flac.FlacFromVTK()
nex = fl.nx - 1;nez = fl.nz - 1

phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
angle = np.zeros(end)
bet = 1.2

rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))

for i in range(137,138):
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)
    trench_ind = np.argmin(z[:,0]) 
    x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]

    m=[]; m2=[]
    x_ocean = mx[phase==phase_oceanic + phase_ecolgite]
    z_ocean = mz[phase==phase_oceanic + phase_ecolgite]
    fig, (bbb,aaa,ccc)= plt.subplots(3,1,figsize=(9,12))    
    start = math.floor(x_trench)
    final = math.floor(np.max(x_ocean))
    bbb.scatter(mx[phase==phase_oceanic + phase_ecolgite],mz[phase==phase_oceanic+phase_ecolgite],color='green')
    bbb.grid()
    bbb.set_ylim(-100,0)
    bbb.set_xlim(start,final)
    bbb.set_aspect('equal')
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    for yy,xx in enumerate(x_grid):
        oz[yy] = np.average(z_ocean[(x_ocean>=px) *(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px) *(x_ocean<=xx)])
        px = xx
    bbb.scatter(ox,oz,color='k')
    kkx=(f2.moving_window_smooth(ox,5))[1:-10]
    kkz=(f2.moving_window_smooth(oz,5))[1:-10]
    kkz=(f2.moving_window_smooth(kkz,5))[1:]
    kkx=kkx[1:]

      
    for kk in range(1,len(kkx)):
        cx1=kkx[kk-1];cx2=kkx[kk]
        cz1=kkz[kk-1];cz2=kkz[kk]
        if (cx2-cx1) != 0:
          m.append((cz2-cz1)/(cx2-cx1))
    qq = kkx[1:]
    for ww in range(1,len(m)):
        cz1=m[ww-1];cz2=m[ww]
        cx1=qq[ww-1];cx1=qq[ww]
        if (cx2-cx1) != 0:
            m2.append((cz2-cz1)/(cx2-cx1))
        else: w1 = ww
    qq2=qq[1:ww]
    mmm=f2.moving_window_smooth(m,5)
    mmm2=f2.moving_window_smooth(m2,6)
    
    aaa.plot(qq,m,color='gray',zorder=1)
    aaa.plot(qq,mmm,color='k',zorder=1)
    ccc.plot([start,final],[0,0],'--',zorder=0,color='red')
    bbb.set_title('frame='+str(i))
    aaa.set_xlim(start,final)
    aaa.grid();ccc.grid() 
    aaa.tick_params(axis='x', labelsize=16)
    aaa.tick_params(axis='y', labelsize=16)
    bbb.tick_params(axis='x', labelsize=16)
    bbb.tick_params(axis='y', labelsize=16)
    ccc.plot(qq2,m2,color='gray')
    ccc.plot(qq2,mmm2,color='k')
    ccc.tick_params(axis='x', labelsize=16)
    ccc.tick_params(axis='y', labelsize=16)
    ccc.set_xlim(start,final)
    ccc.set_ylim(-0.0005,0.0005)
