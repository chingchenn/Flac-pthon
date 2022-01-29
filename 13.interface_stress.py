#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:54:38 2021

@author: jiching
"""
import flac
import os,sys
import numpy as np
import matplotlib.pyplot as plt

model = sys.argv[1]
path = 'D:/model/'+model+'/'
path = '/home/jiching/geoflac/'+model+'/'
#path = '/scratch/jiching/summer2021/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
force=np.zeros(end)
forces=np.zeros(end)
cmap = plt.cm.get_cmap('rainbow')
fig, (ax) = plt.subplots(1,1,figsize=(10,12))
fig3, (ax3) = plt.subplots(1,1,figsize=(10,12))
for flame in range(1,end):
    strain = fl.read_strain(flame)
    sxx = fl.read_sxx(flame)
    strain_rate = fl.read_srII(flame)
    normal_force=fl.read_pres(flame)
    x,z = fl.read_mesh(flame)
    xx = fl.nx
    zz = fl.nz
    element_x = 1/4*(x[:xx-1,:zz-1]+x[1:xx,1:zz]+x[1:xx,:zz-1]+x[:xx-1,1:zz])
    element_z = 1/4*(z[:xx-1,:zz-1]+z[1:xx,1:zz]+z[1:xx,:zz-1]+z[:xx-1,1:zz])
    t = np.zeros(element_x[:,0].shape)
    t[:] = flame*0.2
    ff=np.zeros(nex)
    crust_x = np.zeros(nex)
    crust_z = np.zeros(nex)
    for j in range(nex):
        ind_oceanic = (strain_rate[j,:]>-13.6) 
        if True in ind_oceanic:
            for uu in range(nez):
                height = z[j,uu+1]-z[j,uu]
                forces[flame]+=height*sxx[j,uu]/100
        for mm in range(nez+1):
            ff[j]+=sxx[j,mm-1]*(z[j,mm]-z[j,mm-1])/100
    for kk in range(nez+1):
        maxsxx_id=np.argmin(strain_rate[:,kk-1])
        hight=(z[maxsxx_id,kk]-z[maxsxx_id,kk-1])
        force[flame]+= sxx[maxsxx_id,kk-1]*hight*100
    cb_plot =ax.scatter(element_x[:,0],t,c=ff,cmap=cmap,vmin=1,vmax=100)
    ca_plot = ax3.scatter(element_x[:,0],t,c=strain_rate[:,18],cmap=cmap,vmin=-16,vmax=-12)
ax_cbin = fig.add_axes([0.63, 0.2, 0.23, 0.02])
ax3_cbin = fig3.add_axes([0.63, 0.2, 0.23, 0.02])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ca = fig3.colorbar(ca_plot,cax=ax3_cbin,orientation='horizontal')
ax_cbin.set_title('Sxx (MPa)')
ax.set_xlim(0,1200)
ax.set_ylim(0,fl.time[-1])
ax3_cbin.set_title('Log scale Strain Rate (s^-1)')
ax3.set_xlim(0,1200)
ax3.set_ylim(0,fl.time[-1])
fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxx.jpg')
fig3.savefig('/home/jiching/geoflac/figure/'+model+'_strain_rate.jpg')
# plt.scatter(element_x,element_z,c=sxx,s=2.5)
fig2, (ax2) = plt.subplots(1,1,figsize=(10,6))
# ax2.plot(fl.time,force,c='k')
ax2.plot(fl.time,forces,c='g')
ax2.set_xlim(0,fl.time[-1])
ax2.grid()
fig2.savefig('/home/jiching/geoflac/figure/'+model+'_interface_force.jpg')
