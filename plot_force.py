#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 02:00:10 2022

@author: chingchen
"""

import sys, os
import numpy as np
import flac
import math
import matplotlib
from heapq import nlargest,nsmallest
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.interpolate import  UnivariateSpline,Akima1DInterpolator, PchipInterpolator
#------------------------------------------------------------------------------
# plt.rcParams["font.family"] = "Times New Roman"
# plt.rcParams["figure.figsize"] = (10,12)
#model = sys.argv[1]
model = 'Ref_Cocos'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
path = 'D:/model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T

newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 

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


# def temp_elements(temp):
#     ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
#     return ttt

# def dynamics_pressure(frame):
#     pre = -fl.read_pres(frame) *1e8
#     ooone = pre.flatten()
#     x,z = fl.read_mesh(frame)
#     ele_x,ele_z = flac.elem_coord(x, z)
#     a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
#     fit=(ele_z.flatten()-b)/a
#     dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
#     return x,z,dypre

# def find_slab_median_index2(i):    
#     bet = 1
#     x, z = fl.read_mesh(i)
#     mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
#     ## In this code, we considered the marker phase, not the element phase
#     x_trench = trench_x[i]
#     x_ocean = mx[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]
#     z_ocean = mz[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]
#     # if z_trench> -2 or min(z_ocean)>-200:
#         # continue
#     start = math.floor(x_trench-50)
#     final = math.floor(np.max(x_ocean))
#     x_grid = np.arange(start,final,bet)
#     ox = np.zeros(len(x_grid))
#     oz = np.zeros(len(x_grid))
#     px = start-bet
#     #find initial basalt depth to remove the weage basalt
#     if len(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])==0:
#             print('no')
#     kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
#     x_ocean = x_ocean[z_ocean<kk]
#     z_ocean = z_ocean[z_ocean<kk]
#     # interplate to the grid length "bet"
#     for yy,xx in enumerate(x_grid):
#         if len(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])==0:
#             continue

#         oz[yy] = np.min(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
#         ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
#         px = xx
#     oxx=ox[ox>start]
#     oz=oz[ox>start]
#     ox=oxx
#     return ox,oz

# g = 10
# frame = 140
# dis_range = 25
# depth1 = -200
# depth2 = -10

# bwith = 3
# fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
# Fsbx = np.zeros(4)
# for kkk,tttemp in enumerate([600,700,800,900]):
# # for kkk,tttemp in enumerate([600]):
#     fsb=np.zeros(end)
#     for frame in range(5,end):


#         ###---------------------------------------------------------------------------------------------------------
#         x, z = fl.read_mesh(frame)
#         ele_x, ele_z = flac.elem_coord(x, z)
#         phase = fl.read_phase(frame)
#         _,_,dpre = dynamics_pressure(frame) # N/m^2
#         slab_x,slab_z = find_slab_median_index2(frame)
#         density = fl.read_density(frame)
        
#         xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
#         zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
        
#         z5 =  np.polyfit(xslab,zslab,5)
#         p5 = np.poly1d(z5)
#         w5 = p5(xslab) # 5st poly
#         ###----------------------------------------------------------------------------
#         # colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
#         #           "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
#         #           "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
#         #           "#F67280","#00FF00","#FFFF00","#7158FF"]
#         # phase15= matplotlib.colors.ListedColormap(colors)
#         # cm = plt.cm.get_cmap('seismic_r')
#         # fig, (ax)= plt.subplots(1,1,figsize=(10,6))
#         # # ax.scatter(ele_x,ele_z,c=-dpre/1e6,cmap=cm,vmin=-200, vmax=200,s=40)
#         # ax.pcolormesh(x,-z,density,cmap=cm,vmin=2900,vmax=3500)
#         # # ax.scatter(ele_x,ele_z,c=phase,cmap=phase15,vmin=1, vmax=20,s=100)
#         # ax.set_ylim(300,-0)
#         # ax.set_xlim(200,700)
#         # # ax.set_xlim(trench_x[frame]-200,min(trench_x[frame]+800,1200))
#         # ax.set_aspect('equal')
#         # # xt,zt = fl.read_mesh(frame)
#         # temp = fl.read_temperature(frame)
#         # ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
#         # #plt.scatter(xslab,zslab,c = 'k',s = 50)
#         # # ax.scatter(xslab,w5,c = 'k',s = 30)
#         ###----------------------- Slab sinking force with time-------------------------------
#         #    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 
#         # def slab_sinking_torque(frame):
#             # ------ read data from model -----
#         x, z = fl.read_mesh(frame)
#         ele_x, ele_z = flac.elem_coord(x, z)
#         phase = fl.read_phase(frame)
#         density = fl.read_density(frame)
#         area = fl.read_area(frame)
#         temp = fl.read_temperature(frame)
#         temp_ele = temp_elements(temp)
#         ind_trench = int(trench_index[frame])
#         moment_point_x,moment_point_z  = trench_x[frame], trench_z[frame]
#         # ax.scatter(moment_point_x,-moment_point_z,c = 'green')
#         ref_den = density[-4,:]
#         # ----- empty array and data -----
#         Fsb = 0 
        
        
#         for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
#         # for ii,x_ind in enumerate(range(166,182)):
#             # Choose the eclogite area
#             ind_eclogite = (ele_z[x_ind,:]<-20)*((phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic))
#             if not True in ind_eclogite:
#                 # print(frame,x_ind)
#                 continue
#             top_slab_index = np.where(ind_eclogite)[0][0]
#             litho800 = (temp_ele[x_ind,:]<tttemp)*(ele_z[x_ind,:]<ele_z[x_ind,top_slab_index])
#             if True in litho800:
#                 # ax.scatter(ele_x[x_ind,:][litho800],-ele_z[x_ind,:][litho800],c = 'yellow',s=4)
#                 for ele_index in range(len(ele_x[x_ind,:][litho800])):
#                     x1 = ele_x[x_ind,:][litho800][ele_index]
#                     z1 = ele_z[x_ind,:][litho800][ele_index]
#                     rho_diff = density[x_ind,:][litho800][ele_index]-ref_den[ele_index]
#                     torque_length= (x1-moment_point_x) *1e3        
#                     volume = area[x_ind,:][litho800][ele_index]
#                     Fsb+= torque_length*rho_diff*g*volume
#         fsb[frame]=Fsb
        
        
    
#     ax2.plot(fl.time,fsb,c=newcolors[kkk],label=str(tttemp),lw=4) 
#     ax2.set_xlabel('Time (Myr)',fontsize=16)
#     ax2.set_ylabel('Torque (N)',fontsize=16)
#     ax2.set_xlim(0, fl.time[-1])
#     ax2.tick_params(axis='x', labelsize=16)
#     ax2.tick_params(axis='y', labelsize=16)
#     ax2.grid()
#     ax2.legend(fontsize=20)
#     ax2.spines['bottom'].set_linewidth(bwith)
#     ax2.spines['top'].set_linewidth(bwith)
#     ax2.spines['right'].set_linewidth(bwith)
#     ax2.spines['left'].set_linewidth(bwith)
#     #ax.set_yscale('log')
#     ax2.set_title('Forces of '+model,fontsize=20)
def moving_window_smooth(A,window_width):
    MM = np.zeros(len(A))    
    for kk in range(0,len(A)):
        if kk>(len(A)-window_width):
            MM[kk] = A[kk]
        else:
            MM[kk] = sum(A[kk:kk+window_width])/window_width
    return MM
fig9, (ax)= plt.subplots(1,1,figsize=(10,6))  
bwith=3

## PLOT Torque 
time,fsb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
sb = moving_window_smooth(fsb,3)/1e18
tt = moving_window_smooth(fsu,3)/1e18
ax.plot(time,sb,c='#c06c84',label='slab pull (N)',lw=3)
ax.plot(time,tt,c="#355c7d",label='suction force (N)',lw=3)

#================================figure setting================================
ax.set_ylabel('Torque (10$^{18}$ N)',fontsize=20)
# ax.legend(fontsize=16,loc='upper left')
ax.set_xlabel('Time (Myr)',fontsize=20)

# ax.legend(fontsize=20,facecolor='white',loc='upper left')
ax.grid()
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_xlim(0,30)
ax.set_ylim()
        
    