#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:09:22 2023

@author: chingchen
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
fontsize=30
labelsize=25
bwith = 5
fig6,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(18,10))
axdep1 = ax.twinx()
depth_limit = 160

pres = np.linspace(0,7,100)
sss=np.zeros(len(pres))
TTT = 800-3.5e-8*(pres*3e4-62)**2
lab1=ax.plot(TTT,pres,c='#D14309',lw=5,label='chlorite-peridotite')
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
#ax.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.95, 0.15))    

axdep2 = ax2.twinx()
ax2.vlines(x=650, ymin=0, ymax=3, colors='#6B8E23', ls='-', lw=5,)
ax2.hlines(y=3, xmin=200, xmax=650, colors='#6B8E23', ls='-', lw=5,)

depth = np.linspace(0,300000,100)
sss=np.zeros(len(depth))
for q,dd in enumerate(depth):
    ss1 = 680+0.6e-3*(dd-140e3)
    ss2 = 930-313*(1-np.exp(-dd/7e3))
    sss[q] = max(ss1,ss2)
lab3=axdep2.plot(sss,depth/1e3,c='#FF9900',lw=5,label='solidus')

axdep3 = ax3.twinx()
tempr = np.linspace(0,513)
y = -0.0375 * tempr + 20.1
basalt_change = (50/255, 200/255, 180/255)
ax3.plot(tempr,y,c=basalt_change,lw=5)

x = np.linspace(515,1000)
y = 0.0022 * x - 0.3
lab3=ax3.plot(x,y,c=basalt_change,lw=5,label='basalt-eclogite')
## Dash line for eclogite 
x2 = np.linspace(1000,1200)
y2 = 0.0022 * x2 - 0.3
ax3.plot(x2,y2,c=basalt_change,lw=5,linestyle=(5, (3, 1)))

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
#ax3.legend(lns, labs, fontsize = fontsize-7, loc='lower right',bbox_to_anchor=(0.9, 0.4))

#ax2.text(770,1.5,'schist',fontsize=36)
#ax3.text(570,3.5,'Eclogite',fontsize=26)
ax.set_ylabel('Pressure (GPa)',fontsize=fontsize)
axdep3.set_ylabel('Depth (km)',fontsize=fontsize)
xmajor_ticks=np.array([300,600,900,1200])
ymajor_ticks=np.array([40,80,120,160,160])
for aa in [ax,ax2,ax3]:
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in')
    aa.set_xlim(300,1200)
    aa.grid()
    aa.set_ylim(depth_limit*1000*10*3300/1e9,0)
    aa.set_xlabel('Temperature ($^\circ$C)',fontsize=fontsize)
    aa.set_xticks(xmajor_ticks)
for axdep in [axdep1,axdep2,axdep3]:
    axdep.set_ylim(depth_limit,0)
    axdep.tick_params(axis='y',labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in')
    axdep.set_ylim(depth_limit,0)
    axdep.set_yticks(ymajor_ticks)
    #axdep.grid()
fig6.savefig('/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/phase_diagram_v2.pdf')        