#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:04:18 2022

@author: chingchen
"""
import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
figpath = '/Users/chingchen/Desktop/figure/'
figpath ='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/ThesisNTU/figures/'
bwith = 3
#-------------------------------- setting type --------------------------------
plot_strength                 = 1
save_strength                 = 0 
plot_temperature_profile      = 0
save_temperature_profile      = 0
adiabatic_add                 = 0


#----------------------------------- FIGURE -----------------------------------
if plot_strength:
    #------------------------------- setting figure -------------------------------
    withregion = 0
    strength_fill = 0
    layerz = (0, 12e3, 30e3, 35e3)
    phase=[2,6,4,4]
    tem=3
    #---------------------- define strain rate & Temperature ----------------------
    edot = 1e-14  # high strain rate
    edot = 1e-15  # low strain rate
    deepz = layerz[-1] * 20
    z = np.linspace(0, deepz, num=50000)
    T = f2.continental_geothermal_T3(z,20,6,40)
    
    #---------------------------- read phase from csv -----------------------------
    pu=[]
    for yy in range(20):
        pu.append(f2.phase_pro(yy))
    #----------------------------- creat Dfc array --------------------------------
    pp=[]
    dfc=[0,10,12]
    for qqq in phase:
        for nnn in dfc:    
            pp.append(pu[qqq][nnn])
    Dfc=np.array(pp).reshape(len(phase),3)
    #----------------------------- creat nAE array --------------------------------
    pp=[]
    nae=[3,4,5]
    for qqq in phase:
        for nnn in nae:    
            pp.append(pu[qqq][nnn])
    nAEs=np.array(pp).reshape(len(phase),3)
    #------------------------------------------------------------------------------
    # equation soluiton of plastic stress and viscosity
    frico_strength = f2.plastic_stress(z,layerz,Dfc)
    visc = f2.visc_profile(z, T, edot, layerz, nAEs)
    visco_strength = visc* edot *2 #Pa
    
    ###============================== Plot ========================================
    fig2, (ax) = plt.subplots(1,1,figsize=(6,10))
    applied_strength = np.amin((visco_strength,frico_strength),axis=0)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    mm1,=ax.plot(visco_strength/1e6,z/1000,color='r',linestyle='dashed',alpha=0.8,label = 'elasto-viscous',lw=3)
    mm2,=ax.plot(frico_strength/1e6,z/1000,color='b',linestyle='dashed',alpha=0.8,label = 'elasto-plastic',lw=3)
    mm3,=ax.plot(applied_strength/1e6,z/1000,color='k',lw=5,label = 'final stress')
    mm=[mm3,mm2,mm1]
    ax.tick_params(axis='x', labelsize=26)
    ax.tick_params(axis='y', labelsize=26)
    ax.set_title('Rock Strength',fontsize=30)
    ax.set_xlabel('Strength (MPa)',fontsize=26)
    ax.set_ylabel('Depth (km)',fontsize=26)
    ax.set_ylim(100,0)                                     
    ax.set_xlim(0,1500)
    ax.grid()
    ax.legend(mm, [curve.get_label() for curve in mm],fontsize=20)
    
    ## ------------------------------  Elastic Plot  ------------------------------
    if withregion:
        elastic =z[applied_strength==frico_strength]/1000
        for rr in range(1,len(elastic)):
            if elastic[rr] > -100 and abs(elastic[rr-1]-elastic[rr]) <1 :
                ax.axhspan(elastic[rr-1],elastic[rr],facecolor='royalblue', alpha=0.45)
    if save_strength:
        fig2.savefig(figpath+'strength.pdf')
        
if plot_temperature_profile:
    layerz = (0, 12e3, 30e3, 35e3)
    deepz = layerz[-1] * 20
    z = np.linspace(0, deepz, num=50000)
    TCO = f2.half_space_cooling_T(z, 10, 1330, 40)
    TMO = f2.half_space_cooling_T(z, 10, 1330, 15)
    TC = f2.continental_geothermal_T4(z, 10,1330, 130)
    TM = f2.continental_geothermal_T3(z,20,6,40)
    adiabatic = 3e-5*(1330-10)*10 * adiabatic_add
    CtempC = z/1000*adiabatic+TC
    CtempM = z/1000*adiabatic+TM
    OtempC = z/1000*adiabatic+TCO
    OtempM = z/1000*adiabatic+TMO
    
    ###============================== Plot ========================================
    fig, (ax) = plt.subplots(1,2,figsize=(15,12))
    
    ax[0].plot(OtempC,z/1000,color='#4682B4',label='Oceanic geothermal',lw=5)
    ax[1].plot(OtempM,z/1000,color='#4682B4',label='Oceanic geothermal',lw=5)
    ax[0].plot(CtempC,z/1000,color='#B22222',label='Continental geothermal',lw=5)
    ax[1].plot(CtempM,z/1000,color='#B22222',label='Continental geothermal',lw=5)
    ax[0].set_ylabel('Depth (km)',fontsize=26)
    ax[0].set_title('   Geothermal Gradient \n of Nazca model',fontsize=30)
    ax[1].set_title('   Geothermal Gradient \n of Cocos model',fontsize=30)
    
    for qq in range(len(ax)):
        ax[qq].set_xlim(0,2000)
        ax[qq].set_ylim(150,0)
        # ax3.axes.yaxis.set_visible(False)
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=26)
        ax[qq].tick_params(axis='y', labelsize=26)
        depthmajor_ticks = np.linspace(0,150,num=11)
        ax[qq].set_yticks(depthmajor_ticks)
        ax[qq].set_xlabel('Temperature ($^\circ$C)',fontsize=26)
        ax[qq].grid()
        ax[qq].legend(fontsize=20)
    if save_temperature_profile:
        fig.savefig(figpath+'temperature_profile.pdf')
    