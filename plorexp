#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:50:47 2021

@author: ji-chingchen
"""
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
magma=80
dt=3
lam_tdep = np.logspace(-3,-5,20)
fig, (ax)= plt.subplots(1,1,figsize=(10,8))

def lam(lam0,lam_tdep,delT=1200-10):
    lam = lam0*(np.exp(-lam_tdep*delT))
    return lam
def magma_all(dt,magma,lam,delT=1200-10):
    total_magma= magma * (1-lam * dt)
    delta_fmagma = magma * lam * dt
    return total_magma,delta_fmagma

color = cm.get_cmap('rainbow',20)
color_list = color(np.linspace(0,1,20))
range_lam0 = np.logspace(-11,-12,20)
for kk,lam0 in enumerate (range_lam0):
    for lam_ttdep in lam_tdep:
        llt = lam(lam0,lam_ttdep)
        ax.scatter(lam_ttdep,llt,c=color_list[kk],label=str(range_lam0[kk]))
        # ax.scatter(lam0,llt,c=color_list[kk],label=str(lam0))
ax.set_xscale("log")
# ax.set_yscale("log")
ax.set_xlabel("lambda freeze tdep",fontsize=20)
# ax.set_xlabel("lambda freeze 0",fontsize=20)
ax.set_ylabel("lambda",fontsize=20)
# ax.legend() 

fig2,ax2= plt.subplots(1,1,figsize=(10,8))
for kk,lam0 in enumerate (range_lam0):
    for lam_ttdep in lam_tdep:
        lam_set = lam(lam0,lam_ttdep)
        total_magma,delta_fmagma = magma_all(dt,magma,lam_set)
        ax2.scatter(lam_ttdep,delta_fmagma,c=color_list[kk],label=str(range_lam0[kk]))
        # ax2.scatter(lam0,llam,c=color_list[kk],label=str(lam0))
        print(total_magma)
# ax2.set_xscale("log")
ax2.set_yscale("log")

ax2.set_xlabel("lambda freeze tdep",fontsize=20)
# ax2.set_xlabel("lambda freeze 0",fontsize=20)
ax2.set_ylabel("total magma",fontsize=20)
ax.tick_params(axis='x', labelsize=16 )
ax2.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16 )
ax2.tick_params(axis='y', labelsize=16 )
