#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 14:52:19 2022

@author: ji-chingchen
"""

import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
fontsize=30
labelsize=25
bwith = 5
fig1 = 0 # basalt and eclogit 
fig2 = 0 # perdotite and serpentinite
fig3 = 1 # sediment to schist
if fig1:
###===============================basalt and eclogit ==========================
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    axdep = ax.twinx()
    depth_limit = 200
    depth = np.linspace(0,500e3,1000)
    T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='15 Ma Oceanic geothermal')
    # T = f2.half_space_cooling_T(depth, 10, 1330, 40)+0.4*depth/1e3
    # axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='OC 40 Ma')
    # T = f2.continental_geothermal_T3(depth,20,6,40)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Conten')
    T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')



    tempr = np.linspace(0,513)
    y = -0.0375 * tempr + 20.1
    basalt_change = (50/255, 200/255, 180/255)
    ax.plot(tempr,y,c=basalt_change,lw=8)
    
    x = np.linspace(515,1000)
    y = 0.0022 * x - 0.3
    lab3=ax.plot(x,y,c=basalt_change,lw=8,label='basalt-eclogite')

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
    lab4=ax.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')
    ax.set_ylim(depth_limit*1000*10*3300/1e9,0)
    ax.set_xlim(200,1600)
    
    
    
    lns = lab1+lab2+lab3+lab4#+lab5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(1.7, 0.7))
    
    axdep.set_ylim(depth_limit,0)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep.set_ylabel('Depth (km)',fontsize=fontsize)
    ax.text(240,2.5,'Basalt',fontsize=26)
    ax.text(470,3.5,'Eclogite',fontsize=26)
    # ax.text(480,2.5,'Phase boundary from',fontsize=fontsize)
    # ax.text(520,2.1,'Hacker et al. (2003)',fontsize=fontsize)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    
    
    # ax.set_title('metamorphosed MORB ',fontsize=25)
    fig.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'basalt_phase_diagram'+'.pdf')
    # fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/'+'basalt_phase_diagram'+'.pdf')
###==================== perdotite and serpentinite =============================
if fig2:
    fig2,ax = plt.subplots(1,1,figsize=(10,10))
    depth_limit = 200
    axdep = ax.twinx()
    depth = np.linspace(0,500e3,1000)
    T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='15 Ma Oceanic geothermal')
    # T = f2.half_space_cooling_T(depth, 10, 1330, 40)+0.4*depth/1e3
    # axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='OC 40 Ma')
    # T = f2.continental_geothermal_T3(depth,20,6,40)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Conten')
    T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')
    
    #### KATZ 2003, dry solidus
    x=np.array([0,1,2.6,4,6,7])
    y=np.array([1180,1200,1460,1560,1660,1750])
    pres = np.linspace(0,2000,1000)


    a,b,c = np.polyfit(y, x, 2)
    yyyy = a*pres**2+b*pres+c
    # ax.scatter(y,x,c='k',s=100)
    lab3=ax.plot(pres[pres>1000],yyyy[pres>1000],c='green',lw=5,label='dry solidus')
    # ax.set_xlim(1000,2000)
    # ax.set_ylim(9,0)
    ####
    
    phase_change = (140/255, 20/255, 70/255)
    # fig2,ax = plt.subplots(1,1,figsize=(10,14))
    x = np.linspace(500,730)
    tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
    lab4=ax.plot(x,tt1,c=phase_change,label = 'perdotite-serpentinite',lw=5)
    x = np.linspace(620,730)
    ttold = 2.1 + (0.2-2.1) * (x-730)/(650-730)
    ax.plot(x,ttold,c=phase_change,lw=5)
    
    
    solidus = np.zeros(len(depth))
    for i,kk in enumerate(depth):
        if kk > 113e3:
            solidus[i] = 800+(kk-114)/1e3*5e-3
        else:
            solidus[i] = 800+3.14e-8*(kk-113e3)**2
    
    lab5=axdep.plot(solidus,depth/1e3,c='#FF9900',lw=4,label = 'solidus')
    
    
    lns = lab1+lab2+lab3+lab4+lab5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(1.7, 0.65))
    
    ax.set_ylim(depth_limit*1000*10*3300/1e9,0)
    ax.set_xlim(200,1600)
    axdep.set_ylim(depth_limit,0)
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
    # ax.text(80,3.5,'serpentinite',fontsize=36)
    # ax.text(870,6.5,'perdotite',fontsize=36)
    # ax.legend(fontsize = fontsize-7, loc='upper left')
    fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'serpentinite_phase_diagram'+'.pdf')
    
    ###====================sediment to schist =============================
if fig3:    
    fig3,ax = plt.subplots(1,1,figsize=(10,10))
    axdep = ax.twinx()
    depth_limit = 200
    depth = np.linspace(0,500e3,1000)
    T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='15 Ma Oceanic geothermal')
    # T = f2.half_space_cooling_T(depth, 10, 1330, 40)+0.4*depth/1e3
    # axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='OC 40 Ma')
    # T = f2.continental_geothermal_T3(depth,20,6,40)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Conten')
    T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')

    
    
    
    depth = np.linspace(0,300000,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(depth):
        # ss1 = 620+0.4e-3*(dd-140e3)
        # ss2 = 930-313*(1-np.exp(-dd/18e3))
        ss1 = 680+0.6e-3*(dd-140e3)
        ss2 = 930-313*(1-np.exp(-dd/7e3))
        sss[q] = max(ss1,ss2)
    lab3=axdep.plot(sss,depth/1e3,c='#FF9900',lw=5,label='solidus')
    
    depth = np.linspace(0,300000,100)
    ssss=np.zeros(len(depth))
    for qq,ddd in enumerate(depth):
        ss1 = 620+0.4e-3*(ddd-140e3)
        ss2 = 930-313*(1-np.exp(-ddd/18e3))
        # ss1 = 680+0.6e-3*(ddd-140e3)
        # ss2 = 930-313*(1-np.exp(-ddd/7e3))
        ssss[qq] = max(ss1,ss2)
    lab4=axdep.plot(ssss,depth/1e3,c='green',lw = 5,linestyle = '--',label='sediment-schist')
    
    lns = lab1+lab2+lab3#+lab4#+lab5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(1.7, 0.75))
    
    
    
    ax.set_ylim(depth_limit*1000*10*3300/1e9,0)
    axdep.set_ylim(depth_limit,0)
    ax.set_xlim(200,1600)
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
    # ax.text(80,3.5,'sediment',fontsize=36)
    ax.text(770,1.5,'schist',fontsize=36)
    fig3.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'sediment_phase_diagram'+'.pdf')
    #/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures