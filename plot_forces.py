#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:18:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import matplotlib
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["figure.figsize"] = (10,12)

path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/03model/'
#path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'
###----------------------- Slab sinking force with time-------------------------------

model1 = 'Nazca_aa06'
model2 = 'Nazca_ab04'
model3 = 'Nazca_ab06'
model4 = 'Nazca_ab05'
label1 = '1 P0'
label2 = '10 km ridge'
label3 = '3580'
label4 = '12.5 km ridge'
bwith = 3
fontsize=20
save = 0
color1 = '#858465'
color2 = '#607c85'
color3 = '#871c0f'
color4 = '#303a54'
color_gravity = '#c06c84'
color_suction = '#355c7d'


fig8=1 # suction compare 
fig9=0 # gravity compare 
fig10=0
fig8_save=0
fig9_save=0
def moving_window_smooth(A,window_width):
    MM = np.zeros(len(A))    
    for kk in range(0,len(A)):
        if kk>(len(A)-window_width):
            MM[kk] = A[kk]
        else:
            MM[kk] = sum(A[kk:kk+window_width])/window_width
    return MM

if fig8:
    fig8, (mm2)= plt.subplots(1,1,figsize=(12,5))
    ## PLOT model1 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model1+'_forces.txt').T
    tt = moving_window_smooth(fsu,1)/1e19
    mm2.plot(time,tt,c=color1,label=label1,lw=4)
    
    ## PLOT model2 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model2+'_forces.txt').T
    tt = moving_window_smooth(fsu,1)/1e19
    mm2.plot(time,tt,c=color2,label=label2,lw=4)
    
    ## PLOT model3 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model3+'_forces.txt').T
    tt = moving_window_smooth(fsu,1)/1e19
    mm2.plot(time,tt,c=color3,label=label3,lw=4)
    
    ## PLOT model4 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model4+'_forces.txt').T
    tt = moving_window_smooth(fsu,1)/1e19
    mm2.plot(time,tt,c=color4,label=label4,lw=4)
    #================================figure setting================================
    for aaa in [mm2]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
        aaa.set_xlim(0,40)
    
    mm2.set_ylim(-1,3)
    mm2.set_ylabel('suction torque (10$^{19}$ N$\cdot$m/m)',fontsize=fontsize-4)
    mm2.legend(fontsize=fontsize-6,loc='upper left')
    mm2.set_xlabel('time (Myr)',fontsize=fontsize)
    if fig8_save:
        fig8.savefig(figpath+'fig6b.pdf')
    
    
if fig9:
    fig9, (mm2)= plt.subplots(1,1,figsize=(12,5))
    ## PLOT model1 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model1+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    mm2.plot(time,tt,c=color1,label=label1,lw=4)
    

    ## PLOT model2 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model2+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    mm2.plot(time,tt,c=color2,label=label2,lw=4)
    

    ## PLOT model3 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model3+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    mm2.plot(time,tt,c=color3,label=label3,lw=4)
    
    
    ## PLOT model4 Torque 
    time,fsb,fsu = np.loadtxt(savepath+model4+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    mm2.plot(time,tt,c=color4,label=label4,lw=4)
    su = moving_window_smooth(fsu,1)/1e19
    #mm2.plot(time,su,c=color_suction,lw=4)
    #================================figure setting================================
    for aaa in [mm2]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
        aaa.set_xlim(0,40)

    mm2.set_ylim(0,3)
    mm2.set_ylabel('gravity torque (10$^{19}$ N$\cdot$m/m)',fontsize=fontsize-4)
    mm2.legend(fontsize=fontsize-6,loc='upper left')
    mm2.set_xlabel('time (Myr)',fontsize=fontsize)
    
    if fig9_save:
        fig9.savefig(figpath+'gravity_test.pdf')

if fig10:
    fig10, (mm2,mm3)= plt.subplots(2,1,figsize=(12,10))

    ## PLOT model1 Torque 
    model = model2
    time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    su = moving_window_smooth(fsu,1)/1e19
    mm2.plot(time,tt,c=color_gravity,lw=4)
    mm2.plot(time,su,c=color_suction,lw=4)
    mm2.set_title(model,fontsize=fontsize)
    
    
    ## PLOT model2 Torque 
    model = model3
    time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    tt = moving_window_smooth(fsb,1)/1e19
    su = moving_window_smooth(fsu,1)/1e19
    mm3.plot(time,tt,c=color_gravity,label='gravity',lw=4)
    mm3.plot(time,su,c=color_suction,label='suction',lw=4)
    mm3.set_title(model,fontsize=fontsize)

    
    #================================figure setting================================
    for aaa in [mm2,mm3]:
        aaa.tick_params(labelsize=fontsize)
        aaa.grid()
        for axis in ['top','bottom','left','right']:
            aaa.spines[axis].set_linewidth(bwith)
        aaa.set_xlim(0,40)
        aaa.set_ylim(-1,3)
    mm2.set_ylabel( 'torque (10$^{19}$ N$\cdot$m/m)',fontsize=fontsize-4)
    mm3.legend(fontsize=fontsize-6,loc='upper left')
    mm3.set_xlabel('time (Myr)',fontsize=fontsize)
    
