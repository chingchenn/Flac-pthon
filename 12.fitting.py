# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 00:00:03 2021

@author: JiChing Chen
"""


import flac
import os,sys
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as f2

fig = 1
temp2 = np.loadtxt('D:\\OneDrive - 國立台灣大學\\resarch\\slab 2.0/slab/Gline.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
x,y,z = data.T
sx = x[0]
sy = y[0]
new_cord=np.zeros(len(x))
for uu in range(1,len(x)):
    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)

mindepth=-250

x=new_cord[z>mindepth]
z=z[z>mindepth]
new_cord=x
kk=3
ss=0.001
mean_list=np.zeros(len(x))
median_list=np.zeros(len(x))
tck = interpolate.splrep(x,z,k=kk,s=ss)
zz0=interpolate.splev(x,tck,der=0)    
zz1=interpolate.splev(x,tck,der=1)    
zz2=interpolate.splev(x,tck,der=2) 
yders = interpolate.spalde(x, tck)

    
print('k=',kk,'s=',ss,'mean=',np.mean(zz0-z),'median=',np.median(zz0-z))
    
    

z1=np.polyfit(new_cord,z,1)
p1=np.poly1d(z1)
w1=p1(new_cord)
z2=np.polyfit(new_cord,z,2)
p2=np.poly1d(z2)
w2=p2(new_cord)
z3=np.polyfit(new_cord,z,3)
p3=np.poly1d(z3)
w3=p3(new_cord)
z4=np.polyfit(new_cord,z,4)
p4=np.poly1d(z4)
w4=p4(new_cord)

#cal residual
residual1 = np.sum((w1-z)/len(z))
residual2 = np.sum((w2-z)/len(z))
residual3 = np.sum((w3-z)/len(z))
residual4 = np.sum((w4-z)/len(z))
#p4 
z4=np.polyfit(new_cord,z,5)
p4=np.poly1d(z4)
w4=p4(new_cord)

p4d1=np.polyder(p4,1)
p4d2=np.polyder(p4,2)
w4d1=p4d1(new_cord)
w4d2=p4d2(new_cord)
fig2,(q1,q2,q3)= plt.subplots(3,1,figsize=(7,15))
q1.plot(new_cord,w4,c='k',lw=3)
q1.plot(x,zz0,c='r')

q1.scatter(new_cord,z,c='cyan',s=20)
q1.set_aspect('equal', adjustable='box')  
q2.plot(new_cord,w4d1,c='k')
q2.plot(x,zz1,c='r')

q3.plot(new_cord,w4d2,c='k')
q3.plot(x,zz2,c='r')

q3.set_ylim(-0.03,0.01)
q1.set_ylim(-300,0)
q1.grid();q2.grid();q3.grid() 
q1.tick_params(axis='x', labelsize=16)
q1.tick_params(axis='y', labelsize=16)
q2.tick_params(axis='x', labelsize=16)
q2.tick_params(axis='y', labelsize=16)
q3.tick_params(axis='y', labelsize=16)
q3.tick_params(axis='x', labelsize=16)
# fig2.savefig(path+'frame='+str(i)+'_fig2.png')

