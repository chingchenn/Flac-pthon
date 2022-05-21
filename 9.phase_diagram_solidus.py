#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 14:52:19 2022

@author: ji-chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
fontsize=30
labelsize=25
bwith = 5
fig1 = 1 # basalt and eclogit 
fig2 = 0 # perdotite and serpentinite
fig3 = 0 # sediment to schist
if fig1:
###===============================basalt and eclogit ==========================
    x = np.linspace(0,514)
    y = -0.0375 * x + 20.1
    basalt_change = (50/255, 200/255, 180/255)
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.plot(x,y,c=basalt_change,lw=8)
    x = np.linspace(515,1000)
    y = 0.0022 * x - 0.3
    ax.plot(x,y,c=basalt_change,lw=8,label='basalt-eclogite')

    pressure=np.linspace(0,9,100)
    sss=np.zeros(len(pressure))
    for q,dd in enumerate(pressure):
        if dd<1:
            ss=1050-420*(1-np.exp(-dd*3.3))
        elif dd>2.38:
            ss=(dd+14)*43
        else:
            ss=630+26*((-dd)**2)/2 
        sss[q] = ss
    ax.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')
    ax.set_ylim(0,9)
    ax.set_xlim(0,1200)
    axdep = ax.twinx()
    axdep.set_ylim(0,300)
    ax.legend(fontsize=fontsize-7)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep.set_ylabel('Depth (km)',fontsize=fontsize)
    ax.text(80,3.5,'Basalt',fontsize=36)
    ax.text(570,5.5,'Eclogite',fontsize=36)
    # ax.text(480,2.5,'Phase boundary from',fontsize=fontsize)
    # ax.text(520,2.1,'Hacker et al. (2003)',fontsize=fontsize)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    
    
    # ax.set_title('metamorphosed MORB ',fontsize=25)
    # fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/ThesisNTU/figures/'+'basalt_phase_diagram'+'.pdf')
    ###==================== perdotite and serpentinite =============================
if fig2:
    fig2,ax = plt.subplots(1,1,figsize=(10,10))
    phase_change = (140/255, 20/255, 70/255)
    # fig2,ax = plt.subplots(1,1,figsize=(10,14))
    x = np.linspace(500,730)
    tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
    ax.plot(x,tt1,c=phase_change,label = 'perdotite-serpentinite',lw=5)
    x = np.linspace(670,730)
    ttold = 2.1 + (0.6-2.1) * (x-730)/(670-730)
    ax.plot(x,ttold,c=phase_change,lw=5)
    
    depth = np.linspace(0,900000,100)
    solidus = np.zeros(len(depth))
    for i,kk in enumerate(depth):
        if kk > 80e3:
            solidus[i] = 800
        else:
            solidus[i] = 800+6.2e-8*(kk-80e3)**2
    axdep = ax.twinx()
    axdep.plot(solidus,depth/1e3,c='pink',lw=4,label = 'solidus')
    ax.set_ylim(0,9)
    ax.set_xlim(0,1200)
    axdep.set_ylim(0,300)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep.set_ylabel('Depth (km)',fontsize=fontsize)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.text(80,3.5,'serpentinite',fontsize=36)
    ax.text(770,5.5,'perdotite',fontsize=36)
    ax.legend(fontsize = fontsize-7)
    
    ###====================sediment to schist =============================
if fig3:    
    fig3,ax = plt.subplots(1,1,figsize=(10,10))
    depth = np.linspace(0,300000,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(depth):
        ss1 = 680+0.6e-3*(dd-140e3)
        ss2=930-313*(1-np.exp(-dd/7e3))
        sss[q] = max(ss1,ss2)
    axdep = ax.twinx()
    axdep.plot(sss,depth/1e3,c='#FF9900',lw=5)
    ax.set_ylim(0,9)
    axdep.set_ylim(0,300)
    ax.set_xlim(0,1200)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep.set_ylabel('Depth (km)',fontsize=fontsize)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.text(80,3.5,'sediment',fontsize=36)
    ax.text(770,5.5,'schist',fontsize=36)
    
    