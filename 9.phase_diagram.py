#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 09:48:29 2021

@author: ji-chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
x = np.linspace(0,514)
y = -0.0375 * x + 20.1
basalt_change = (50/255, 200/255, 180/255)

fig,ax = plt.subplots(1,1,figsize=(10,10))
ax.plot(x,y,c=basalt_change)
ax.set_ylim(0,5)
ax.set_xlim(0,1000)

x = np.linspace(500,1000)
y = 0.0022 * x - 0.3
ax.plot(x,y,c=basalt_change)
# ax.grid()
ax.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16 )
ax.set_xlabel('Temperature ($^\circ$C)',fontsize=20)
ax.set_ylabel('Pressure (Gpa)',fontsize=20)
 
ax.plot([500,500],[0.8,6],'k--')
# ax.set_title('metamorphosed MORB ',fontsize=25)
# fig.savefig('/Users/ji-chingchen/Desktop/'+'basalt_phase_diagram'+'.pdf')
phase_change = (140/255, 20/255, 70/255)
fig2,ax2 = plt.subplots(1,1,figsize=(10,14))
x = np.linspace(500,730)
tt1 = 2.1 +(7.5-2.1)* (x - 730)/ (500-730)
ax2.plot(x,tt1,c=phase_change)
ax2.set_ylim(0,8)
ax2.set_xlim(450,800)
x = np.linspace(500,730)
tt2 = 2.1 + (0.2-2.1) * (x-730)/(500-730)
ax2.plot(x,tt2,c=phase_change)
x = np.linspace(650,730)
ttold = 2.1 + (0.2-2.1) * (x-730)/(650-730)
ax2.plot(x,ttold,'k--')
ax2.tick_params(axis='x', labelsize=16 )
ax2.tick_params(axis='y', labelsize=16 )
ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=20)
ax2.set_ylabel('Pressure (Gpa)',fontsize=20)

