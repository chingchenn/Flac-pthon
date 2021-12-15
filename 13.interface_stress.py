# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:54:38 2021

@author: jiching
"""
import flac
import os
import numpy as np
import matplotlib.pyplot as plt

model='s1303'
path = 'D:/model/'+model+'/'
# path = '/home/jiching/geoflac/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
force=np.zeros(end+1)
forces=np.zeros(end+1)
cmap = plt.cm.get_cmap('rainbow')
fig, (ax) = plt.subplots(1,1,figsize=(10,12))
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
        ind_oceanic = (strain_rate[j,:]>-12.6) 
        if True in ind_oceanic:
            for uu in range(nez):
                height = z[j,uu+1]-z[j,uu]
                forces[flame]+=height*sxx[j,uu]*1000/10e6
        for mm in range(nez+1):
            ff[j]+=sxx[j,mm-1]*(z[j,mm]-z[j,mm-1])*1000/10e6
    for kk in range(nez+1):
        maxsxx_id=np.argmin(strain_rate[:,kk-1])
        hight=(z[maxsxx_id,kk]-z[maxsxx_id,kk-1])
        force[flame]+= sxx[maxsxx_id,kk-1]*hight*1000/10e6
    cb_plot =ax.scatter(element_x[:,0],t,c=ff,cmap=cmap)
ax_cbin = fig.add_axes([0.93, 0.18, 0.23, 0.02])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ax_cbin.set_title('Stress')
ax.set_xlim(0,1200)
ax.set_ylim(0,50)
# fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxx.jpg')
# plt.scatter(element_x,element_z,c=sxx,s=2.5)
fig2, (ax2) = plt.subplots(1,1,figsize=(10,6))
# ax2.plot(fl.time,force,c='k')
ax2.plot(fl.time,forces,c='g')
ax2.set_xlim(0,50)
# fig2.savefig('/home/jiching/geoflac/figure/'+model+'_interface_force.jpg')