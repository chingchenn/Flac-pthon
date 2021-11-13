#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 17:01:45 2021

@author: jiching
"""

import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"

layerz = (0, 2e3, 6e3, 16e3)   # 1st elem must be 0
Dfc = ((2800, 11,4e7),  #phase  11
        (2880,30,4e7), #phase 3
        (3200,30,4e7), #phase 16
        (3300,30,4e7)) #phase  4
         
nAEs = ((3.00, 5.00e+2, 2.00e+5),
        (3.05, 1.25e-1, 3.76e+5),
        (3.00, 7.00e+4, 5.20e+5),
        (3.00, 7.00e+4, 5.20e+5))
#layerz = (0, 6e3, 16e3)   # 1st elem must be 0
#Dfc = ( (2880,30,4e7), #phase 3
#        (3200,30,4e7), #phase 16
#        (3300,30,4e7)) #phase  4
         
#nAEs = ((3.05, 1.25e-1, 3.76e+5),
#        (3.00, 7.00e+4, 5.20e+5),
#        (3.00, 7.00e+4, 5.20e+5))
# layerz = (0, 18e3, 30e3, 40e3)   # 1st elem must be 0
# Dfc = ((2800,30,4e7), #phase 2
#         (2900,30,4e7), #phase 6
#         (3200,30,4e7), #phase19
#         (3300,30,4e7)) #phase 4
# nAEs = ((3.05, 1.25e-1, 2.76e+5),
#         (3.05, 1.25e-1, 3.76e+5),
#         (3.1, 7.00e+5, 5.76e+5),
#         (3.00, 7.00e+5, 5.76e+5))
#edot = 1e-14  # high strain rate
edot = 1e-15  # low strain rate
deepz = layerz[-1] * 10
z = np.linspace(0, deepz, num=1000)

#------------------------------------------------------------------------------
# equation soluiton of plastic stress and viscosity
frico_strength = f2.plastic_stress(z,layerz,Dfc)
#T = f2.continental_geothermal_T(z,20,6,45)
T = f2.half_space_cooling_T(z, 10, 1330, 15)
visc = f2.visc_profile(z, T, edot, layerz, nAEs)
visco_strength = visc* edot *2 #Pa
#------------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(6,10))
applied_strength = np.amin((visco_strength,frico_strength),axis=0)
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
mm1,=ax.plot(visco_strength/1e6,-z/1000,'--r',alpha=0.8,label = 'visco')
mm2,=ax.plot(frico_strength/1e6,-z/1000,'--b',alpha=0.8,label = 'plastic')
mm3,=ax.plot(applied_strength/1e6,-z/1000,'k',lw=2,label = 'final stress')
ax.set_ylim(-100,0)                                     
ax.set_xlim(0,1e3)
ax.tick_params(axis='x', labelsize=26)
ax.tick_params(axis='y', labelsize=26)
ax.set_title('Rock Strength',fontsize=30)
ax.set_xlabel('Strength (MPa)\nTemperature ($^\circ$C)',fontsize=26)
ax.set_ylabel('Depth (km)',fontsize=26)
ax.grid()
ax2 = ax.twinx()
temp = z/1000*0.6+T
mm4,=ax2.plot(temp,-z/1000,color='green',label='temperature')
mm=[mm1,mm2,mm3,mm4]
ax2.set_xlim(0,1500)
ax2.set_ylim(-100,0)
ax2.axes.yaxis.set_visible(False)
ax.legend(mm, [curve.get_label() for curve in mm],fontsize=20)
ax2.set_xlabel('Temperature ($^C)',fontsize=26)

# ax.fill_between(applied_strength/1e6, -z/1000, facecolor='tab:cyan', interpolate=True,alpha=0.6)
# fig.savefig('/home/jiching/geoflac/figure/'+'strength_profile'+'.png')
## Intergal
tol = 0
tt=np.zeros(len(z))
for ww in range(1,len(z)):
    tol += applied_strength[ww] *(z[ww]-z[ww-1])/1e9
    tt[ww] = applied_strength[ww]*(z[ww]-z[ww-1])/1e9
print(tol)
