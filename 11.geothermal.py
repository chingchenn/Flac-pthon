#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 20:54:20 2021

@author: ji-chingchen
"""

import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt

#plt.rcParams["font.family"] = "Times New Roman"

layerz = (0, 2e3, 6e3, 16e3)   # 1st elem must be 0
deepz = layerz[-1] * 20
z = np.linspace(0, deepz, num=1000)
depth=170
#------------------------------------------------------------------------------
# equation soluiton of plastic stress and viscosity

T_con = f2.continental_geothermal_T3(z,20,6,45)
T_oce = f2.half_space_cooling_T(z, 10, 1330, 40)
T_con4 = (1330-0)/depth*z/1000
#------------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(9,15))
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)

ax.tick_params(axis='x', labelsize=26)
ax.tick_params(axis='y', labelsize=26)
ax.set_title('Geotherm',fontsize=30)
ax.set_xlabel('Temperature ($^\circ$C)',fontsize=26)
ax.set_ylabel('Depth (km)',fontsize=26)
ax.grid()
for ii in range(len(T_con)):
    if T_con4[ii]>=1330:
        T_con4[ii]=1330
    if T_oce[ii]>=1330:
        T_oce[ii]=1330

T_oce = z/1000*0.6+T_oce
T_con = z/1000*0.6+T_con
T_con4 = z/1000*0.6+T_con4
ax.plot(T_con,-z/1000,color='#8B008B',lw=6,label='continental geotherm')
ax.plot(T_oce,-z/1000,color='#4169E1',lw=6,label='oceanic geotherm')
ax.plot(T_con4,-z/1000,color='g',lw=6,label='one layer continental geotherm')
ax.axvline(1330,-z[-1]/1000,0,c='r')
ax.set_xlim(0,2000)
ax.set_ylim(-z[-1]/1000,0)
# ax.legend(fontsize=26)
ax.legend(fontsize=16)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
fig.savefig('/home/jiching/geoflac/figure/goethermal.png')
