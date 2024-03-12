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
fig3 = 0 # sediment to schist
fig4 = 0 # chlorite
fig5 = 1 # chlorite, perdotite and serpentinite
fig6 = 1 # ALL three figures
if fig1:
###===============================basalt and eclogit ==========================
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    axdep = ax.twinx()
    depth_limit = 200
    # depth = np.linspace(0,500e3,1000)
    # T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    # lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='-',label='15 Ma Oceanic geothermal')
    # # T = f2.half_space_cooling_T(depth, 10, 1330, 40)+0.4*depth/1e3
    # # axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='OC 40 Ma')
    # # T = f2.continental_geothermal_T3(depth,20,6,40)+0.4*depth/1e3
    # # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Conten')
    # T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')



    tempr = np.linspace(0,513)
    y = -0.0375 * tempr + 20.1
    basalt_change = (50/255, 200/255, 180/255)
    ax.plot(tempr,y,c=basalt_change,lw=5)
    
    x = np.linspace(515,1000)
    y = 0.0022 * x - 0.3
    lab3=ax.plot(x,y,c=basalt_change,lw=5,label='basalt-eclogite')

    pressure_limit = 2.45 # GPa
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
    ax.set_xlim(200,1200)
    
    # x = np.linspace(710,1050)
    # y = -1.25/350*x+5
    # ax.plot(x,y,c='#FF9900',lw=5) # solidus

    # x = np.linspace(680,1050)
    # y = 0.65/400*x-0.45625
    # ax.plot(x,y,c='#FF9900',lw=5) # solidus
    
    # lns = lab1+lab2+lab3+lab4#+lab5
    lns = lab3+lab4#+lab5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.9, 0.4))
    
    axdep.set_ylim(depth_limit,0)
    axdep.tick_params(axis='y', labelsize=labelsize)
    ax.tick_params(axis='x', labelsize=labelsize)
    ax.tick_params(axis='y', labelsize=labelsize)
    ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep.set_ylabel('Depth (km)',fontsize=fontsize)
    # ax.text(240,2.5,'Basalt',fontsize=26)
    ax.text(570,3.5,'Eclogite',fontsize=26)
    # ax.text(480,2.5,'Phase boundary from',fontsize=fontsize)
    # ax.text(520,2.1,'Hacker et al. (2003)',fontsize=fontsize)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    
    
    # ax.set_title('metamorphosed MORB ',fontsize=25)
    # fig.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'basalt_phase_diagram'+'.pdf')
    # fig.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Method/'+'basalt_phase_diagram.pdf')
    # fig.savefig('/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/'+'basalt_phase_diagram'+'.pdf')
###==================== perdotite and serpentinite =============================
if fig2:
    fig2,ax = plt.subplots(1,1,figsize=(15,10))
    depth_limit = 200
    axdep = ax.twinx()
    # depth = np.linspace(0,500e3,1000)
    # T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    # lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='-',label='15 Ma Oceanic geothermal')
    # # T = f2.half_space_cooling_T(depth, 10, 1330, 40)+0.4*depth/1e3
    # # axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='OC 40 Ma')
    # # T = f2.continental_geothermal_T3(depth,20,6,40)+0.4*depth/1e3
    # # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Conten')
    # T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')
    
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
    lab4=ax.plot(x,tt1,c=phase_change,label = 'serpentinite-perdotite',lw=5)
    x = np.linspace(620,730)
    ttold = 2.1 + (0.2-2.1) * (x-730)/(650-730)
    ax.plot(x,ttold,c=phase_change,lw=5)
    
    
    depth = np.linspace(0,300000,100)
    presss = np.linspace(0,9.9,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(presss):
        ss1 = 980+0.6e-3*(dd*3.3e4-140e3)
        ss2 = 1090-178*(1-np.exp(-dd*4.125))
        sss[q] = max(ss1,ss2)
    lab5=ax.plot(sss,presss,c='#FF9900',lw=5,label='solidus')

    
    
    lns = lab3+lab4+lab5#+lab1+lab2+
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
    # fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'serpentinite_phase_diagram'+'.pdf')
    # fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Method/'+'serpentinite_phase_diagram.pdf')
    
    ###====================sediment to schist =============================
if fig3:    
    fig3,ax = plt.subplots(1,1,figsize=(7,10))
    axdep = ax.twinx()
    depth_limit = 200
    # depth = np.linspace(0,500e3,1000)
    # T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    # lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle=':',label='15 Ma Oceanic geothermal')
    # T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    # lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')

    ax.vlines(x=650, ymin=0, ymax=3, colors='green', ls='-', lw=5,)
    ax.hlines(y=3, xmin=200, xmax=650, colors='green', ls='-', lw=5,)
    
    
    depth = np.linspace(0,300000,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(depth):
        ss1 = 680+0.6e-3*(dd-140e3)
        ss2 = 930-313*(1-np.exp(-dd/7e3))
        sss[q] = max(ss1,ss2)
    lab3=axdep.plot(sss,depth/1e3,c='#FF9900',lw=5,label='solidus')
    
    # depth = np.linspace(0,300000,100)
    # ssss=np.zeros(len(depth))
    # for qq,ddd in enumerate(depth):
    #     ss1 = 620+0.4e-3*(ddd-140e3)
    #     ss2 = 930-313*(1-np.exp(-ddd/18e3))
    #     # ss1 = 680+0.6e-3*(ddd-140e3)
    #     # ss2 = 930-313*(1-np.exp(-ddd/7e3))
    #     ssss[qq] = max(ss1,ss2)
    # lab3=axdep.plot(ssss,depth/1e3,c='green',lw = 5,linestyle = '--',label='sediment-schist')
    
    # lns = lab1+lab2+lab3#+lab4#+lab5
    # labs = [l.get_label() for l in lns]
    # ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(1.7, 0.75))
    
    
    
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
    # fig3.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'sediment_phase_diagram'+'.pdf')

    
        ###====================phase 16 =============================
if fig4:    
    fig4,ax = plt.subplots(1,1,figsize=(10,10))
    axdep = ax.twinx()
    depth_limit = 200
    depth = np.linspace(0,500e3,1000)
    T = f2.half_space_cooling_T(depth, 10, 1330, 15)+0.4*depth/1e3
    lab1=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='-',label='15 Ma Oceanic geothermal')
    T = f2.continental_geothermal_T4(depth, 10,1330, 130)+0.4*depth/1e3
    lab2=axdep.plot(T,depth/1e3,c='gray',lw = 3,linestyle='--',label='Continental geothermal')
    
    pres = np.linspace(0,7,100)
    sss=np.zeros(len(pres))
    # for q,dd in enumerate(pres):
    TTT = 800-3.5e-8*(pres*3e4-62)**2
    lab3=ax.plot(TTT,pres,c='#D14309',lw=5,label='chlorite-peridotite')

    lns = lab1+lab2+lab3
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
    # fig4.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'chlorite_phase_diagram'+'.pdf')
    # fig4.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Method/'+'chlorite_phase_diagram.pdf')
if fig5:
    fig5,ax = plt.subplots(1,1,figsize=(10,10))
    axdep = ax.twinx()
    depth_limit = 200
    
    pres = np.linspace(0,7,100) # GPa
    TTT=np.zeros(len(pres))
    for q,dd in enumerate(pres):
        ss1 = 700+150*dd
        ss2 = 764.545+53.7*dd
        ss3 = 973.447-43.478*dd
        ss4 = -63.64*dd+1040.016
        ss5 = -200*dd+1640
        TTT[q] = min(ss1,ss2,ss3,ss4,ss5)
    ax.scatter(700,0,c='#D14309',s=50)
    ax.scatter(800,0.6,c='#D14309',s=50)
    ax.scatter(880,2.15,c='#D14309',s=50)
    ax.scatter(830,3.3,c='#D14309',s=50)
    ax.scatter(760,4.4,c='#D14309',s=50)
    ax.scatter(600,5.2,c='#D14309',s=50)
    
    # TTT = 800-3.5e-8*(pres*3e4-62)**2
    # TTT = 700-4.3e-8*(pres*1.8e4-62)**2
    lab1=ax.plot(TTT,pres,c='#D14309',lw=5,label='chlorite-peridotite')
    #### KATZ 2003, dry solidus
    x=np.array([0,1,2.6,4,6,7])
    y=np.array([1180,1200,1460,1560,1660,1750])
    pres = np.linspace(0,2000,1000)


    a,b,c = np.polyfit(y, x, 2)
    yyyy = a*pres**2+b*pres+c
    # ax.scatter(y,x,c='k',s=100)
    lab3=ax.plot(pres[pres>1000],yyyy[pres>1000],c='green',lw=5,label='dry solidus')
    
    ####
    
    phase_change = (140/255, 20/255, 70/255)
    # fig2,ax = plt.subplots(1,1,figsize=(10,14))
    x = np.linspace(500,730)
    tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
    lab4=ax.plot(x,tt1,c=phase_change,label = 'serpentinite-perdotite',lw=5)
    x = np.linspace(620,730)
    ttold = 2.1 + (0.2-2.1) * (x-730)/(650-730)
    ax.plot(x,ttold,c=phase_change,lw=5)
    
    
    # depth = np.linspace(0,300000,100)
    presss = np.linspace(0,9.9,100)
    sss=np.zeros(len(presss))
    for q,dd in enumerate(presss):
        ss1 = 980+0.6e-3*(dd*3.3e4-140e3)
        ss2 = 1090-178*(1-np.exp(-dd*4.125))
        sss[q] = max(ss1,ss2)
    lab5=ax.plot(sss,presss,c='#FF9900',lw=5,label='solidus')

    
    
    lns = lab3+lab4+lab5+lab1#+lab2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.95, 0.15))
    
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
    # fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'+'serpentinite_phase_diagram'+'.pdf')
    # fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Method/'+'serpentinite_phase_diagram.pdf')
if fig6: 
    fig6,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(18,10))
    axdep1 = ax.twinx()
    depth_limit = 150
    
    pres = np.linspace(0,7,100)
    sss=np.zeros(len(pres))
    TTT1 = 800-3.5e-8*(pres*3e4-62)**2
    lab1=ax.plot(TTT1,pres,c='#D14309',lw=5,label='chlorite-peridotite')
    #### KATZ 2003, dry solidus
    x=np.array([0,1,2.6,4,6,7])
    y=np.array([1180,1200,1460,1560,1660,1750])
    pres = np.linspace(0,2000,1000)

    a,b,c = np.polyfit(y, x, 2)
    yyyy = a*pres**2+b*pres+c
    lab3=ax.plot(pres[pres>1000],yyyy[pres>1000],c='green',lw=5,label='dry solidus')
    
    ####
    
    phase_change = (140/255, 20/255, 70/255)
    x = np.linspace(500,730)
    tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
    lab4=ax.plot(x,tt1,c=phase_change,label = 'serpentinite-perdotite',lw=5)
    x = np.linspace(620,730)
    ttold = 2.1 + (0.2-2.1) * (x-730)/(650-730)
    ax.plot(x,ttold,c=phase_change,lw=5)
    
    depth = np.linspace(0,300000,100)
    presss = np.linspace(0,9.9,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(presss):
        ss1 = 980+0.6e-3*(dd*3.3e4-140e3)
        ss2 = 1090-178*(1-np.exp(-dd*4.125))
        sss[q] = max(ss1,ss2)
    lab5=ax.plot(sss,presss,c='#FF9900',lw=5,label='solidus')
    lns = lab3+lab4+lab5+lab1#+lab2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.95, 0.15))    

    axdep2 = ax2.twinx()
    ax2.vlines(x=650, ymin=0, ymax=3, colors='green', ls='-', lw=5,)
    ax2.hlines(y=3, xmin=200, xmax=650, colors='green', ls='-', lw=5,)
    
    depth = np.linspace(0,300000,100)
    sss=np.zeros(len(depth))
    for q,dd in enumerate(depth):
        ss1 = 680+0.6e-3*(dd-140e3)
        ss2 = 930-313*(1-np.exp(-dd/7e3))
        sss[q] = max(ss1,ss2)
    lab3=axdep2.plot(sss,depth/1e3,c='#FF9900',lw=5,label='solidus')
    
    # ----------------------------basalt-eclogite--------------------
    axdep3 = ax3.twinx()
    tempr = np.linspace(0,513)
    y = -0.0375 * tempr + 20.1
    basalt_change = (50/255, 200/255, 180/255)
    ax3.plot(tempr,y,c=basalt_change,lw=5)
    
    x = np.linspace(515,1000)
    y = 0.0022 * x - 0.3
    lab3=ax3.plot(x,y,c=basalt_change,lw=5,label='basalt-eclogite')
    
    # ---------------------------------solidus----------------------------------
    pressure_limit = 2.45 # GPa
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
    lab4=ax3.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')

    lns = lab3+lab4
    labs = [l.get_label() for l in lns]
    ax3.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.9, 0.4))

    ax2.text(770,1.5,'schist',fontsize=36)
    ax3.text(570,3.5,'Eclogite',fontsize=26)
    ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
    axdep3.set_ylabel('Depth (km)',fontsize=fontsize)
    xmajor_ticks=np.array([300,600,900,1200])
    for aa in [ax,ax2,ax3]:
        aa.spines['bottom'].set_linewidth(bwith)
        aa.spines['top'].set_linewidth(bwith)
        aa.spines['right'].set_linewidth(bwith)
        aa.spines['left'].set_linewidth(bwith)
        aa.tick_params(axis='x', labelsize=labelsize)
        aa.tick_params(axis='y', labelsize=labelsize)
        aa.set_xlim(250,1250)
        aa.set_ylim(depth_limit*1000*10*3300/1e9,0)
        aa.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
        aa.set_xticks(xmajor_ticks)
    for axdep in [axdep1,axdep2,axdep3]:
        axdep.set_ylim(depth_limit,0)
        axdep.tick_params(axis='y', labelsize=labelsize-5)
        axdep.set_ylim(depth_limit,0)
    # fig6.savefig('/Users/chingchen/Desktop/Eclogite_flat_slab/phase_diagram.pdf')        