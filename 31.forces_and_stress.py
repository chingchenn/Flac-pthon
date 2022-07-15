#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
model = sys.argv[1]
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
time,trench_index, trench_x, trench_z = np.loadtxt('/home/jiching/geoflac/data/trench_for_'+str(model)+'.txt').T
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
for i in range(1,end,20):
    dis, ssxz = shearstress_indistance(i)
    ax2.plot(dis-trench_x[i],ssxz,c = newcolors[i],lw=4)
ax2.set_xlim(0,700)
#ax2.set_xlim(np.average(trench_x),1200)
ax2.set_xlabel('Distance (km)',fontsize=16)
ax2.set_ylabel('$\sigma_{xz}$ (MPa)',fontsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.grid()
bwith = 3
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)

fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxx_time.png')
    

#------------------------------------------------------------------------------
if __name__ == '__main__':
    fsb=np.zeros(end)
    ft=np.zeros(end)
    for i in range(1,end):
        fsb[i] = slab_sinking_force(i)
        ft[i] = mantle_traction_force(i)

    fs.save_3txt(model+'_forces','/home/jiching/geoflac/data/',fl.time,fsb,ft)
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
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    #ax.set_yscale('log')
    ax.set_title('Forces of '+model,fontsize=20)
    fig.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force.png')
