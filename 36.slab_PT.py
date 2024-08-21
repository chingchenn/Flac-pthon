#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 22:40:44 2022

@author: chingchen
"""

import flac
import math
import sys, os
import matplotlib
import numpy as np
from matplotlib import cm
from scipy import interpolate
import function_for_flac as fd
# import Main_creat_database as Md
import function_savedata as fs
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"

Cocos = 1
Nazca = 0
if Cocos:
    max_dis = 620
    xmin = 500
    xmax = 900
    model = 'Ref_Cocos'
if Nazca:
    max_dis = 515
    xmin = 300
    xmax = 1000
    model = 'Nazca_a0702'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
path = '/Users/chingchen/Desktop/model/'
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1

phase_uppercrust = 2
phase_oceanic = 3
phase_mantle1 = 4
phase_schist = 5
phase_mantle2 = 8
phase_serpentinite = 9
phase_sediment = 10
phase_sediment_1 = 11
phase_eclogite = 13
phase_lowercrust = 14
phase_hydratedmantle = 16
phase_oceanic_1 = 17
phase_eclogite_1 = 18


frame = 150
labelsize = 20
bwith = 3
colors = ["#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#F67280","#00FF00","#FFFF00","#7158FF"]

# phase_color= ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
#       "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
#       "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
#       "#F67280","#00FF00","#FFFF00","#7158FF"]
# phase15= matplotlib.colors.ListedColormap(phase_color)

colors2=[
 '#C98F49', '#92C0DF', '#2553C7', '#FFFFFF', '#6495ED',
 '#2E8B57', '#524B52', '#9A32CD', '#6B8E23','#D4DBF4',
 '#D8BFD8','#999999','#F2C85B','#92C0DF','#999999',
 '#4CC552','#999999','#999999','#999999','#999999']
phase8= matplotlib.colors.ListedColormap(colors2)

def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt
# ------ read data from model -----
x, z = fl.read_mesh(frame)
ele_x, ele_z = flac.elem_coord(x, z)
temp = fl.read_temperature(frame)
temp_ele = temp_elements(temp)
phase = fl.read_phase(frame)
pressure = fl.read_pres(frame) * -0.1
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
ind_trench = int(trench_index[frame])

# ----- empty array and data -----
slab_P = np.zeros(len(ele_z))
slab_T = np.zeros(len(ele_z))
slab_z = np.zeros(len(ele_z))
slab_P_bot = np.zeros(len(ele_z))
slab_T_bot = np.zeros(len(ele_z))
slab_z_bot = np.zeros(len(ele_z))
## ---------------------------------------------
fig1,ax = plt.subplots(1,1,figsize=(10,10))
ax.pcolormesh(x,-z,phase,cmap=phase8,vmin=1, vmax=20)
ax.set_xlim(xmin,xmax)
ax.set_ylim(150,0)
ax.set_aspect('equal')
# ----- Start Calculation ------
for ii,x_ind in enumerate(range(ind_trench-10,len(ele_z))):
    ind_slab = (ele_z[x_ind,:]> -300)*(ele_z[x_ind,:]< -5)*((phase[x_ind,:] == phase_eclogite)+(phase[x_ind,:] == phase_oceanic))
    if not True in ind_slab:
        continue
    if ele_x[x_ind,0]<max_dis:
        top_slab_index = np.where(ind_slab)[0][-1]
    elif ele_x[x_ind,0]<max_dis+300:
        top_slab_index = np.where(ind_slab)[0][-1]
    else:
        continue
    slab_T[x_ind] = temp_ele[x_ind,top_slab_index]
    slab_z[x_ind] = -ele_z[x_ind,top_slab_index]
    slab_P[x_ind] = pressure[x_ind,top_slab_index] # GPa
    ax.scatter(ele_x[x_ind,top_slab_index],-ele_z[x_ind,top_slab_index],c='k',s=5)
    


for ii,x_ind in enumerate(range(ind_trench-10,len(ele_z))):
    ind_slab = (ele_z[x_ind,:]> -300)*(ele_z[x_ind,:]< -5)*(phase[x_ind,:] == phase_sediment)
    if not True in ind_slab:
        continue
    if ele_x[x_ind,0]<max_dis:
        bot_slab_index = np.where(ind_slab)[0][-1]
    elif ele_x[x_ind,0]<max_dis+300:
        bot_slab_index = np.where(ind_slab)[0][-1]
    else:
        continue
    
    slab_T_bot[x_ind] = temp_ele[x_ind,bot_slab_index]
    slab_z_bot[x_ind] = -ele_z[x_ind,bot_slab_index]
    slab_P_bot[x_ind] = pressure[x_ind,bot_slab_index] # GPa
    
    ax.scatter(ele_x[x_ind,bot_slab_index],-ele_z[x_ind,bot_slab_index],c='r',s=5)



slab_P2 = np.zeros(len(ele_z))
slab_T2 = np.zeros(len(ele_z))
slab_z2 = np.zeros(len(ele_z))
slab_P2_bot = np.zeros(len(ele_z))
slab_T2_bot = np.zeros(len(ele_z))
slab_z2_bot = np.zeros(len(ele_z))
for ii,z_ind in enumerate(range(0,len(ele_z[0]))):
    ind_slab = (ele_z[:,z_ind]> -300)*(ele_z[:,z_ind]< -5)*((phase[:,z_ind] == phase_eclogite)+(phase[:,z_ind] == phase_oceanic)+(phase[:,z_ind] == phase_sediment))
    if not True in ind_slab:
        continue
    top_slab_index = np.where(ind_slab)[0][-1]
    bot_slab_index = np.where(ind_slab)[0][0]
    
    slab_T2_bot[z_ind] = temp_ele[bot_slab_index,z_ind]
    slab_z2_bot[z_ind] = -ele_z[bot_slab_index,z_ind]
    slab_P2_bot[z_ind] = pressure[bot_slab_index,z_ind] # GPa
    ax.scatter(ele_x[bot_slab_index,z_ind],-ele_z[bot_slab_index,z_ind],c='k',s=5)


## bottom of sediment
for ii,z_ind in enumerate(range(0,len(ele_z[0]))):
    ind_slab = (ele_z[:,z_ind]> -300)*(ele_z[:,z_ind]< -5)*(phase[:,z_ind] == phase_sediment)
    if not True in ind_slab:
        continue
    top_slab_index = np.where(ind_slab)[0][0]

    slab_T2[z_ind] = temp_ele[top_slab_index,z_ind]
    slab_z2[z_ind] = -ele_z[top_slab_index,z_ind]
    slab_P2[z_ind] = pressure[top_slab_index,z_ind] # GPa
    ax.scatter(ele_x[top_slab_index,z_ind],-ele_z[top_slab_index,z_ind],c='r',s=5)
   
           
fig2,ax2 = plt.subplots(1,1,figsize=(10,10))

x = np.linspace(0,514)
y = -0.0375 * x + 20.1
basalt_change = (50/255, 200/255, 180/255)
ax2.plot(x,y,c=basalt_change,lw=8)
x = np.linspace(515,1300)
y = 0.0022 * x - 0.3
ax2.plot(x,y,c=basalt_change,lw=8,label='basalt-eclogite')
pressure_limit = 7 # GPa
pressure=np.linspace(0,pressure_limit,100)


pressure_limit2 = 2.45
pressure=np.linspace(0,pressure_limit2,100)
sss=np.zeros(len(pressure))
for q,dd in enumerate(pressure):
    if dd<1:
        ss=1050-420*(1-np.exp(-dd*3.3))
    elif dd>2.38:
        ss=(dd+14)*43
    else:
        ss=630+26*((-dd)**2)/2 
    sss[q] = ss
lab4=ax2.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')
# ax.set_ylim(depth_limit*1000*10*3300/1e9,0)
# ax.set_xlim(200,1600)


# sss=np.zeros(len(pressure))
# for q,dd in enumerate(pressure):
#     if dd<1:
#         ss=1050-420*(1-np.exp(-dd*3.3))
#     elif dd>2.38:
#         ss=(dd+14)*43
#     else:
#         ss=630+26*((-dd)**2)/2 
#     sss[q] = ss
# ax2.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')
# x = np.linspace(710,1050)
# y = -1.25/350*x+5
# ax2.plot(x,y,c='#FF9900',lw=5) # solidus

# x = np.linspace(680,1050)
# y = 0.65/400*x-0.45625
# ax2.plot(x,y,c='#FF9900',lw=5) # solidus
ax2.set_ylim(0,pressure_limit)
ax2.set_ylim(0,8)
ax2.set_xlim(0,1400)
axdep = ax2.twinx()
axdep.set_ylim(0,pressure_limit*1e9/3300/10/1e3)
# ax2.legend(fontsize=labelsize,loc='upper left')
axdep.tick_params(axis='y', labelsize=labelsize)
ax2.tick_params(axis='x', labelsize=labelsize)
ax2.tick_params(axis='y', labelsize=labelsize)
ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=labelsize)
ax2.set_ylabel('Pressure (GPa)',fontsize=labelsize)
axdep.set_ylabel('Depth (km)',fontsize=labelsize)
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)

axdep.scatter(slab_T[slab_z>0],slab_z[slab_z>0],c=colors[10],s=10)
axdep.scatter(slab_T2[slab_z2>0],slab_z2[slab_z2>0],c='purple',s=10)
axdep.scatter(slab_T2_bot[slab_z2_bot>0],slab_z2_bot[slab_z2_bot>0],c=colors[10],s=10)
axdep.scatter(slab_T_bot[slab_z_bot>0],slab_z_bot[slab_z_bot>0],c='purple')
#ax2.scatter(slab_T[slab_P>0],slab_P[slab_P>0],c='#455E45')

# fig2.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Discussion/'+'Mexico_slab_PT.pdf')
