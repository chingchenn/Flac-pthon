#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:07:54 2024

@author: chingchen
"""

import numpy as np
import function_for_flac as f2
import matplotlib.pyplot as plt


datapath = '/Users/chingchen/Desktop/data/'
# temp2=np.loadtxt(datapath+'g4_table4.txt')
# data = temp2[~np.isnan(temp2).any(axis=1)]
# x,y,z = data.T
# sx = x[0]
# sy = y[0]
# mindepth=-200
# new_cord=np.zeros(len(x))
# for uu in range(1,len(x)):
#     new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
# x=new_cord[z>mindepth]
# z=z[z>mindepth]
# x=x[z<-10]
# z=z[z<-10]
# new_cord=x

# fig0,(q1)= plt.subplots(1,1,figsize=(10,8))
# q1.plot(new_cord,z,'k--')
# q1.set_aspect('equal', adjustable='box')
# q1.set_ylim(mindepth,0)
# q1.grid()
# q1.tick_params(axis='both', labelsize=16)


model = 'b5'
temp2=np.loadtxt(datapath+model+'_table4.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
x,y,z = data.T
sx = x[0]
sy = y[0]
mindepth=-350
new_cord=np.zeros(len(x))
for uu in range(1,len(x)):
    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
x=new_cord[z>mindepth]
z=z[z>mindepth]
x=x[z<-10]
z=z[z<-10]
new_cord=x

fig1,(q2)= plt.subplots(1,1,figsize=(18,15))
q2.set_title(model,fontsize = 20)
q2.plot(new_cord,z,'k',lw=3)
q2.set_xlim(0,900)
q2.set_aspect('equal', adjustable='box')
q2.set_ylim(mindepth,0)
# q2.grid()
q2.tick_params(axis='both', labelsize=16)
xmajor_ticks = np.linspace(900,0,num=10)
q2.set_xticks(xmajor_ticks)

# fig2,(q3)= plt.subplots(1,1,figsize=(18,15))
# temp2=np.loadtxt(datapath+model+'_table4.txt')
# data = temp2[~np.isnan(temp2).any(axis=1)]
# x,y,z = data.T
# sx = x[0]
# sy = y[0]
# mindepth=-350

# q3.set_title(model,fontsize = 20)
# q3.plot(x,z,'k',lw=3)

# new_cord=np.zeros(len(x))
# for uu in range(1,len(x)):
#     new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
# x=new_cord[z>mindepth]
# z=z[z>mindepth]
# x=x[z<-10]
# z=z[z<-10]
# new_cord=x
# q4 = q3.twiny()
# q4.plot(new_cord,z,'k',lw=3)
# q4.set_aspect('equal')
# q3.set_ylim(mindepth,0)
# q4.set_xlim(0,1200)
# q3.tick_params(axis='both', labelsize=16)
# # xmajor_ticks = np.linspace(900,0,num=10)
# # q2.set_xticks(xmajor_ticks)