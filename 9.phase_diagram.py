#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 09:48:29 2021

@author: ji-chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
fig1 = 0
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
    ax.set_ylim(0,8)
    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=26)
    ax.set_ylabel('Pressure (GPa)',fontsize=26)
    ax.text(80,3.5,'Basalt',fontsize=36)
    ax.text(570,5.5,'Eclogite',fontsize=36)
    ax.text(480,2.5,'Phase boundary from',fontsize=26)
    ax.text(520,2.1,'Hacker et al. (2003)',fontsize=26)
    bwith = 5
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
     
    # ax.set_title('metamorphosed MORB ',fontsize=25)
    # fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/ThesisNTU/figures/'+'basalt_phase_diagram'+'.pdf')
    ###====================perdotite and serpentinite =============================
    phase_change = (140/255, 20/255, 70/255)
    # fig2,ax = plt.subplots(1,1,figsize=(10,14))
    x = np.linspace(500,730)
    tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
    ax.plot(x,tt1,c=phase_change,label = 'perdotite-serpentinite')
    x = np.linspace(670,730)
    ttold = 2.1 + (0.6-2.1) * (x-730)/(670-730)
    ax.plot(x,ttold,c=phase_change)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=20)
    ax.set_ylabel('Pressure (Gpa)',fontsize=20)
    
    ###====================sediment to schist =============================
    
    # fig3,ax = plt.subplots(1,1,figsize=(10,14))
    
    depth = np.linspace(0,100000)
    pressure=3000*10*depth/1e9  # (GPa) static pressutr = density(mantle) * g * depth
    dd=20*1e3
    dpressure = 3000*10*dd/1e9
    ax.plot([650,650],[dpressure,8],'r--',label='sediment-schist')
    ax.set_ylim(0,8)
    ax.set_xlim(0,1000)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=20)
    ax.set_ylabel('Pressure (Gpa)',fontsize=20)
    ax.legend(fontsize=20)
##=========================plot seeting ==================================
# fig2,ax2 = plt.subplots(1,1,figsize=(10,10))
# depth = np.linspace(0,120000,100)
# pressure=3000*10*depth/1e9  # (GPa)

# sss=np.zeros(len(depth))
# for q,dd in enumerate(depth):
#     ss1 = 680+0.6e-3*(dd-140e3)
#     ss2=930-313*(1-np.exp(-dd/7e3))
#     sss[q] = max(ss1,ss2)
# ax2.plot(sss,depth/1e3,c='#FF9900',lw=5)
# BB = 24*depth/1e3-190
# with_plot = (BB>600)*(BB<780)
# ax2.plot(BB[with_plot],depth[with_plot]/1e3,c='#FF6699',lw=5)
# AA=-17*depth/1e3+2030
# with_plot = (depth/1e3<82)*(depth/1e3>70)
# ax2.plot(AA[with_plot],depth[with_plot]/1e3,c='#FF6699',lw=5)

# #--------------------------------------------------------------------
# depth = np.linspace(20e3,100000,100)
# rC = np.zeros(len(depth))
# for q,dd in enumerate(depth):
#     rC1=2030/3+7/3*dd/1e3
#     rC2=880
#     rC[q] = min(rC1,rC2)
# # with_plot = (rC>710)*(rC<840)
# ax2.plot(rC,depth/1e3,c='purple',lw=5)


# ax2.set_xlim(0,1500)
# ax2.set_ylim(120,0)
# ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=30)
# ax2.set_ylabel('Depth (km)',fontsize=30)
# bwith = 5
# ax2.spines['bottom'].set_linewidth(bwith)
# ax2.spines['top'].set_linewidth(bwith)
# ax2.spines['right'].set_linewidth(bwith)
# ax2.spines['left'].set_linewidth(bwith)
# ax2.tick_params(axis='x', labelsize=26)
# ax2.tick_params(axis='y', labelsize=26)
# fig2.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/master03/Seminar/my present/sediment_phase_diagram.pdf')

##=========================plot seeting ==================================
fig3,ax3 = plt.subplots(1,1,figsize=(10,10))
depth = np.linspace(0,200000,100)
pressure=3000*10*depth/1e9  # (GPa)

sss=np.zeros(len(pressure))
sss2=np.zeros(len(pressure))
sss3=np.zeros(len(pressure))
ss=0
ss2=0
ss3=0
for q,dd in enumerate(pressure):
    if dd<1:
        ss=1050-420*(1-np.exp(-dd*3.3))
    elif dd>2.7:
        ss=(dd+14)*43
    else:
        ss=630+26*((-dd)**2)/2 
    sss[q] = ss
    sss2[q] = ss2
    sss3[q] = ss3
ax3.plot(sss,pressure,c='#FF9900',lw=5)
x = np.linspace(720,1040)
y = 1.4/(-320) * x + 5.85
ax3.plot(x,y,c='#708090',lw=5)
x = np.linspace(670,1040)
y = 0.6/370 * x + 0.7-670*0.6/370
ax3.plot(x,y,c='#708090',lw=5)
# ax3.plot(sss3,pressure,c='#A80359',lw=5)
# ax3.scatter(700,2.25,c='b',s=300)
# BB = 24*depth/1e3-190
# with_plot = (BB>600)*(BB<780)
# ax3.plot(BB[with_plot],depth[with_plot]/1e3,c='#FF6699',lw=5)
# AA=-17*depth/1e3+2030
# with_plot = (depth/1e3<82)*(depth/1e3>70)
# ax3.plot(AA[with_plot],depth[with_plot]/1e3,c='#FF6699',lw=5)

# #--------------------------------------------------------------------
# depth = np.linspace(20e3,100000,100)
# rC = np.zeros(len(depth))
# for q,dd in enumerate(depth):
#     rC1=2030/3+7/3*dd/1e3
#     rC2=880
#     rC[q] = min(rC1,rC2)
# # with_plot = (rC>710)*(rC<840)
# ax3.plot(rC,depth/1e3,c='purple',lw=5)


ax3.set_xlim(0,1200)
ax3.set_ylim(0,5)
ax3.set_xlabel('Temperature ($^\circ$C)',fontsize=30)
ax3.set_ylabel('Depth (km)',fontsize=30)
bwith = 5
ax3.spines['bottom'].set_linewidth(bwith)
ax3.spines['top'].set_linewidth(bwith)
ax3.spines['right'].set_linewidth(bwith)
ax3.spines['left'].set_linewidth(bwith)
ax3.tick_params(axis='x', labelsize=26)
ax3.tick_params(axis='y', labelsize=26)
# fig3.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/master03/Seminar/my present/sediment_phase_diagram.pdf')