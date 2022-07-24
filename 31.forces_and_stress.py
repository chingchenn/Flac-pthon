#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import matplotlib
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
#model = sys.argv[1]
model = 'b0702k'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Users/ji-chingchen/Desktop/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
###----------------------- Slab sinking force with time-------------------------------
#    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 
def slab_sinking_force(frame):
    phase_mantle1 = 4
    phase_mantle2 = 8
    phase_serpentinite = 9
    phase_hydratedmantle = 16
    phase_eclogite = 13
    phase_eclogite_1 = 18
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    density = fl.read_density(frame)
    area = fl.read_area(frame)
    rho_diff = np.zeros(nez)
    slab_area = np.zeros(nez)
    Fsb = 0; g = 10
    for z_ind in range(1,len(ele_z[0])):
    #for z_ind in range(11,15):
        ind_eclogite = (phase[:,z_ind]==phase_eclogite) + (phase[:,z_ind] == phase_eclogite_1) + (phase[:,z_ind] == phase_hydratedmantle) + (phase[:,z_ind] == phase_mantle2)
        man_eclogite = (phase[:,z_ind]== phase_mantle1) + (phase[:,z_ind] == phase_serpentinite)
        if True in ind_eclogite and True in man_eclogite:
            den_mantle = np.average(density[man_eclogite,z_ind])
            den_eco = np.average(density[ind_eclogite,z_ind])
            rho_diff[z_ind] = den_eco - den_mantle
            slab_area[z_ind] = area[ind_eclogite,z_ind].sum()
        Fsb += rho_diff[z_ind] * g * slab_area[z_ind]
    return Fsb # N/m (2D)

###---------------------- Mantle flow traction force with time -------------------------------
#   Ft = sum(shear stress * dx) = \sigma_xz * dx 
def mantle_traction_force(frame):
    phase_uppercrust = 2
    phase_lowercrust = 14
    x,z, = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x,z)
    sxz = fl.read_sxz(frame)
    phase = fl.read_phase(frame)
    dx = np.zeros(nex)
    stressxz = np.zeros(nex)
    Ft = 0
    ## Base of the upper plate
    for qq in range(1,nex):
        upper_plate = (phase[qq,:]== phase_uppercrust) + (phase[qq,:] == phase_lowercrust)
        if True in upper_plate:
            last_deep= np.argmin(ele_z[qq,upper_plate])
            #print(qq,last_deep)
            dx[qq] = (x[qq+1,last_deep]-x[qq,last_deep])*1e3 # km --> m
            stressxz[qq] = abs(sxz[qq,last_deep]*1e8) # stress kb --> N/m^2
        Ft += stressxz[qq] * dx[qq] # N/m
    return Ft # N/m (2D)
###---------------------- Mantle suction force with time -------------------------------
#def suction_force(frame):
phase_mantle1 = 4
phase_mantle2 = 8
phase_serpentinite = 9
phase_hydratedmantle = 16
phase_eclogite = 13
phase_eclogite_1 = 18
phase_oceanic = 3
phase_oceanic_1 = 17
phase_sediment = 10
phase_sediment_1 = 11
phase_schist = 5
frame = 55

x, z = fl.read_mesh(frame)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(frame)
density = fl.read_density(frame)
pressure = fl.read_pres(frame)
area = fl.read_area(frame)
pdiff = np.zeros(nez)
slab_area = np.zeros(nez)
Psub = 0; Ptop = 0
aatop = 0;aasub = 0
onepre=-pressure.flatten()
a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))*100  # MPa

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
plt.scatter(ele_x,ele_z,c = phase,cmap = phase15,vmin = 1,vmax = 20)
#plt.ylim(-250,-50)
#plt.xlim(600,800)
ins = 15 # pressure area
for x_ind in range(int(trench_index[frame]),len(ele_z)):
#for x_ind in range(145,146):
    #man_eclogite = (phase[x_ind,:]== phase_mantle1) + (phase[x_ind,:] == phase_serpentinite) + (phase[x_ind,:] == phase_hydratedmantle)
    ind_oceanic = (phase[x_ind,:] == phase_oceanic) + (phase[x_ind,:] == phase_eclogite)+(phase[x_ind,:] == phase_oceanic_1) + (phase[x_ind,:] == phase_eclogite_1)
    subducted_sed = (phase[x_ind,:] == phase_sediment) + (phase[x_ind,:] ==phase_sediment_1) + (phase[x_ind,:] == phase_schist)
    ### NEED TO DELETED THE LINEAR PRESSURE HERE. S
    if True in (ind_oceanic+subducted_sed):
        oceanic_plate_index = [ww for ww, x in enumerate(ind_oceanic+subducted_sed) if x]
        bottom_slab_ind = max(oceanic_plate_index) + 1
        top_slab_ind = min(oceanic_plate_index) - 1
        if top_slab_ind  < 0:
            continue
        submantle = (phase[x_ind,bottom_slab_ind:bottom_slab_ind+ins]==phase_mantle1) + (phase[x_ind,bottom_slab_ind:bottom_slab_ind+ins]==phase_mantle2) + (phase[x_ind,bottom_slab_ind:bottom_slab_ind+ins]==phase_serpentinite)
        topmantle = (phase[x_ind,top_slab_ind:top_slab_ind+ins]==phase_mantle1) + (phase[x_ind,top_slab_ind:top_slab_ind+ins]==phase_mantle2) + (phase[x_ind,top_slab_ind:top_slab_ind+ins]==phase_serpentinite)
    
    if True in submantle:
        Psub += dypre[x_ind,bottom_slab_ind:bottom_slab_ind+ins]
        Ptop += dypre[x_ind,top_slab_ind:top_slab_ind+ins]
        aasub += (area[x_ind,bottom_slab_ind:bottom_slab_ind+ins]).sum()
        aatop += (area[x_ind,top_slab_ind:top_slab_ind+ins]).sum()
    print(x_ind,aatop,aasub)
        #plt.scatter(ele_x[x_ind,bottom_slab_ind:bottom_slab_ind+ins],ele_z[x_ind,bottom_slab_ind:bottom_slab_ind+ins], c= 'w')
        #plt.scatter(ele_x[x_ind,top_slab_ind:top_slab_ind+ins],ele_z[x_ind,top_slab_ind:top_slab_ind+ins], c= 'r')     
        #print(x_ind,Fsub-Ftop)   
        
        #rho_diff[z_ind] = den_eco - den_mantle
        #slab_area[z_ind] = area[ind_eclogite,z_ind].sum()
    #Fsb += rho_diff[z_ind] * g * slab_area[z_ind]


def shearstress_indistance(frame):
    phase_uppercrust = 2
    phase_lowercrust = 14
    x,z, = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x,z)
    phase = fl.read_phase(frame)
    sxz = fl.read_sxz(frame)
    stressxz = np.zeros(nex)
    for qq in range(1,nex):
        upper_plate = (phase[qq,:]== phase_uppercrust) + (phase[qq,:] == phase_lowercrust)
        if True in upper_plate:
            last_deep= np.argmin(ele_z[qq,upper_plate])
            stressxz[qq] = sxz[qq,last_deep]*1e2
    ssxz = fd.moving_window_smooth(stressxz,8)
    return ele_x[:,0], ssxz

fig, (ax2) = plt.subplots(1,1,figsize=(15,9))
rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))
for i in range(40,end,20):
    dis, ssxz = shearstress_indistance(i)
    ax2.plot(dis-trench_x[i],ssxz,c = newcolors[i],lw=4,label=str(round(fl.time[i],1))+' Myr')
ax2.set_xlim(0,700)
#ax2.set_xlim(np.average(trench_x),1200)
ax2.set_title(model+' $\sigma_{xz}$ and distance',fontsize = 24)
ax2.set_xlabel('Distance from trench (km)',fontsize=16)
ax2.set_ylabel('$\sigma_{xz}$ (MPa)',fontsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
ax2.legend(fontsize=20)
bwith = 3
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)

#fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxz_dis.png')

#------------------------------------------------------------------------------
if __name__ == '__main__':
    fsb=np.zeros(end)
    ft=np.zeros(end)
    ratio=np.zeros(end)
    for i in range(1,end):
        fsb[i] = slab_sinking_force(i)
        ft[i] = mantle_traction_force(i)
        if ft[i] ==0:
            ratio[i] = 0
        else:
            ratio[i] = fsb[i]/ft[i]

    fs.save_4txt(model+'_forces','/home/jiching/geoflac/data/',fl.time,fsb,ft,ratio)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    sb = fd.moving_window_smooth(fsb[fsb>0],8)
    tt = fd.moving_window_smooth(ft[ft>0],8)
    ax.plot(fl.time[fsb>0],sb,c='#c06c84',label='slab pull (N/m)',lw=4)
    ax.plot(fl.time[ft>0],tt,c="#355c7d",label='traction force (N/m)',lw=4)
    #ax.scatter(fl.time[fsb>0],fsb[fsb>0],c='#c06c84',label='slab pull (N/m)')
    #ax.scatter(fl.time[ft>0],ft[ft>0],c="#355c7d",label='traction force (N/m)')
    ax.legend(fontsize=16,loc='upper left')
    #================================figure setting================================
    ax.set_xlabel('Time (Myr)',fontsize=16)
    ax.set_ylabel('Force (N/m)',fontsize=16)
    ax.set_xlim(0, fl.time[-1])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    #ax.set_yscale('log')
    ax.set_title('Forces of '+model,fontsize=20)
    fig.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force.png')
    fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
    ratio_f = fd.moving_window_smooth(ratio[ratio>0],5)
    ax2.plot(fl.time[ratio>0],ratio_f,c="#355c7d",label='ratio of these forces)',lw=4)
    ax2.set_xlabel('Time (Myr)',fontsize=16)
    ax2.set_xlim(0, fl.time[-1])
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig2.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force_ratio.png')
