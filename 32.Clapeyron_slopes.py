# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:25:05 2022

@author: grace
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt



fig, (ax)= plt.subplots(1,1,figsize=(6,6))
axdep = ax.twinx()

rho = 3450
g = 10

temp = np.linspace(1000,2400,1000)
zzz = (-660e3+(temp-1330)* 500/6) # m
pres = -zzz* rho * g/1e6  #MPa
ax.plot(temp,pres/1e3,lw = 7,c = 'darkgreen') # GPa

labelsize=15;fontsize = 20
axdep.set_ylim(200,800)
ax.set_ylim(200*1e3*rho*g/1e9,800*1e3*rho*g/1e9)
bwith = 3
ax.set_xlim(1000,2400)

ax.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
axdep.set_ylabel('Depth (km)',fontsize=fontsize)
for qq in [ax,axdep]:
    qq.spines['bottom'].set_linewidth(bwith)
    qq.spines['top'].set_linewidth(bwith)
    qq.spines['right'].set_linewidth(bwith)
    qq.spines['left'].set_linewidth(bwith)
    qq.tick_params(axis='x', labelsize=labelsize)
    qq.tick_params(axis='y', labelsize=labelsize)
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(temp,pres)
print(slope)