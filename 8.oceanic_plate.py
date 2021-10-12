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
path = '/Volumes/My Book/model/'+model+'/'
# path = '/Volumes/SSD500/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
# flk = flac.FlacFromVTK()
nex = fl.nx - 1;nez = fl.nz - 1

phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
angle = np.zeros(end)
bet = 2

rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))

for i in range(137,138):
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)
    trench_ind = np.argmin(z[:,0]) 
    x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]

    m=[]; m2=[]
    x_ocean = mx[(phase==phase_ecolgite)+(phase==phase_oceanic)]
    z_ocean = mz[(phase==phase_ecolgite)+(phase==phase_oceanic)]
    fig, (bbb,aaa,ccc)= plt.subplots(3,1,figsize=(9,12))    
    start = math.floor(x_trench)
    final = math.floor(np.max(x_ocean))
    # bbb.scatter(x_ocean,z_ocean,color='green',s=2)
    bbb.grid()
    bbb.set_ylim(-100,0)
    bbb.set_xlim(start,final)
    bbb.set_aspect('equal')
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    #find initial basalt depth to remove the weage basalt
    kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
    x_ocean = x_ocean[z_ocean<kk]
    z_ocean = z_ocean[z_ocean<kk]
    bbb.scatter(x_ocean,z_ocean,color='orange',s=20)
    
    for yy,xx in enumerate(x_grid):
        oz[yy] = np.average(z_ocean[(x_ocean>=px) *(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px) *(x_ocean<=xx)])
        px = xx
    bbb.scatter(ox,oz,color='cyan',s=10)
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
        else: www = ww
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
    # ccc.tick_params(axis='x', labelsize=16)
    ccc.tick_params(axis='y', labelsize=16)
    ccc.set_xlim(start,final)
    ccc.set_ylim(-0.0005,0.0005)


    fig2,(q1,q2,q3)= plt.subplots(3,1,figsize=(9,12))
    ox = ox[oz>-100]
    oz = oz[oz>-100]
    z1=np.polyfit(ox,oz,4)
    p4=np.poly1d(z1)
    w1=p4(ox)
    q1.plot(ox,w1,c='k',lw=3)
    q1.scatter(ox,oz,c='cyan',s=20)
    q1.set_xlim(start,final)
    
    p3=np.polyder(p4,1)
    p2=np.polyder(p4,2)
    
    w2=p3(ox)
    q2.plot(ox,w2,c='k')
    w3=p2(ox)
    q3.plot(ox,w3,c='k')
    q1.set_ylim(-100,0)
    q3.plot([start,final],[0,0],'--',zorder=0,color='red')
    q1.set_title('frame='+str(i))
    q2.set_xlim(start,final)
    q1.grid();q2.grid();q3.grid() 
    q1.tick_params(axis='x', labelsize=16)
    q1.tick_params(axis='y', labelsize=16)
    q2.tick_params(axis='x', labelsize=16)
    q2.tick_params(axis='y', labelsize=16)
    
    # q3.tick_params(axis='x', labelsize=16)
    q3.tick_params(axis='y', labelsize=16)
    q3.set_xlim(start,final)
    # q3.set_ylim(-0.0005,0.0005)
    fig.savefig(path+model+'frame='+str(i)+'_fig1.png')
    fig2.savefig(path+model+'frame='+str(i)+'_fig2.png')