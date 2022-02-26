#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 09:48:29 2021

@author: ji-chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
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
ax.text(480,2.5,'Phase boundary from)',fontsize=26)
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
pressure=3000*10*depth/1e9  # (GPa)
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
