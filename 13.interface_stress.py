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
path = '/home/jiching/geoflac/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
force=np.zeros(end)
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
    for hh in range(nex):
        for mm in range(nez+1):
            ff[hh]+=sxx[hh,mm-1]*(z[hh,mm]-z[hh,mm-1])*1000/10e6
    for kk in range(nez+1):
        maxsxx_id=np.argmin(sxx[:,kk-1])
        hight=z[maxsxx_id,kk]-z[maxsxx_id,kk-1]
        force[flame]+= sxx[maxsxx_id,kk-1]*hight*1000/10e6
    cb_plot =ax.scatter(element_x[:,0],t,c=ff,cmap=cmap)
ax_cbin = fig.add_axes([0.93, 0.18, 0.23, 0.02])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ax_cbin.set_title('Stress')
ax.set_xlim(0,1200)
ax.set_ylim(0,50)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxx.jpg')
# plt.scatter(element_x,element_z,c=sxx,s=2.5)
fig2, (ax2) = plt.subplots(1,1,figsize=(10,6))
ax2.plot(fl.time,force,)
fig2.savefig('/home/jiching/geoflac/figure/'+model+'_interface_force.jpg')