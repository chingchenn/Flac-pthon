#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:35:49 2022

@author: chingchen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
DIR='/Users/chingchen/Desktop/GMT'
ff=pd.read_csv(DIR+'/'+'2020Hu.csv')
list_profile = [0,  # ALU
                82, # RYU
                105,# SAM
                118,# SOC
                ]

name = ['ALU','RYU','SAM','SOC']
xlime_list = [700,400,700,400]
# for kk in range(0,len(ff)):
for uu,kk in enumerate(list_profile):
    index,lat,lon,trench,az=ff.loc[kk].tolist()
    #-------------------------------call gmt to cut the trace----------------------------------
    if index <=9:
        qq = '00'+str(index)
    elif index >= 100:
        qq = str(index)
    else:
        qq = '0'+str(index)
    temp=np.loadtxt(DIR+'/slab_cut_153/'+qq+'_slab_distance.txt')
    data = temp[~np.isnan(temp).any(axis=1)]
    xslab,zslab = data.T
    plt.rcParams["font.family"] = "Times New Roman"
    fig, (aa) = plt.subplots(1,1,figsize=(12,10))
    aa.plot(xslab,-zslab,'k--',lw = 3)
    aa.set_aspect('equal', adjustable='box')
    aa.tick_params(axis='x', direction="in",labelsize=30,pad=20)
    aa.tick_params(axis='y', direction="in",labelsize=30,pad=10)
    aa.tick_params(width=3,length=10)
    aa.set_ylim(200,0)
    aa.set_xlim(0,xlime_list[uu])
    depthmajor_ticks = np.linspace(0,200,num=5)

    
    if xlime_list[uu]==400:
        x_ticks = np.linspace(0,400,num=9)
        bwith = 5
        aa.tick_params(axis='x', direction="in",labelsize=50,pad=20)
        aa.tick_params(axis='y', direction="in",labelsize=50,pad=10)
        aa.set_xticks(x_ticks)
        aa.set_yticks(depthmajor_ticks)
        aa.spines['bottom'].set_linewidth(bwith)
        aa.spines['top'].set_linewidth(bwith)
        aa.spines['right'].set_linewidth(bwith)
        aa.spines['left'].set_linewidth(bwith)
        # fig.savefig(DIR+'/'+name[uu]+'profile_for_slab2.0.pdf')
    else:
        x_ticks = np.linspace(0,700,num=15)
        bwith = 3
        aa.tick_params(axis='x', direction="in",labelsize=30,pad=20)
        aa.tick_params(axis='y', direction="in",labelsize=30,pad=10)
        aa.set_xticks(x_ticks)
        aa.set_yticks(depthmajor_ticks)
        aa.spines['bottom'].set_linewidth(bwith)
        aa.spines['top'].set_linewidth(bwith)
        aa.spines['right'].set_linewidth(bwith)
        aa.spines['left'].set_linewidth(bwith)
        # fig.savefig(DIR+'/'+name[uu]+'profile_for_slab2.0.pdf')
   