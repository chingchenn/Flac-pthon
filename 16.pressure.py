#!/usr/bin/env python

import flac
import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import function_for_flac as f2
import matplotlib.pyplot as plt

#model = sys.argv[1]
model = 'Nazca_a0614'
# frame = int(sys.argv[2])
frame = 78
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
g = 10
cc = plt.cm.get_cmap('RdYlBu_r')
x,z=fl.read_mesh(frame)
phase = fl.read_phase(frame)

def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    new_pre = np.zeros((len(x)-1, len(x[0])-1))
    for xx in range(len(ele_x)):
        den = fl.read_density(frame)
        ref_den = den[-4,:]
        new_pre[xx,:] = pre[xx,:]+ref_den*g*ele_z[xx,:]*1e3
    return x,z,new_pre
x,z,dypre = dynamics_pressure(frame)
fig2, (ax2,ax,ax3)= plt.subplots(3,1,figsize=(36,10))
cb_plot=ax2.pcolormesh(x,z,dypre/1e6,cmap=cc,vmin=-200, vmax=200)
ax2.set_xlim(200,1000)
ax2.set_ylim(-200,0)
ax_cbin = fig2.add_axes([0.43, 0.02, 0.23, 0.03])
cb2 = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
tick_font_size = 20
cb2.ax.tick_params(labelsize=tick_font_size)
ax_cbin.set_title('MPa',fontsize=20)
ax2.set_aspect('equal')



x,z=fl.read_mesh(frame)
phase = fl.read_phase(frame)
ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
pre=fl.read_pres(frame)*1e8 
onepre=-pre.flatten()

pre = -fl.read_pres(frame) *1e8
x,z = fl.read_mesh(frame)
ele_x,ele_z = nodes_to_elements(x, z)
new_pre = np.zeros((len(x)-1, len(x[0])-1))
def dynamics_pressure2(frame):
    pre = -fl.read_pres(frame) *1e8
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    new_pre = np.zeros((len(x)-1, len(x[0])-1))
    for zz in range(len(ele_x[0])):
        den = fl.read_density(frame)
        ref_den = den[:,zz]
        new_pre[:,zz] = pre[:,zz]+ref_den*g*ele_z[:,zz]*1e3
    return x,z,new_pre
# return x,z,new_pre

a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))/1e6 #MPa
x,z,dypre = dynamics_pressure2(frame)
cb_plot=ax.pcolormesh(ele_x,ele_z,dypre/1e6,cmap=cc,vmin=-200, vmax=200)#,vmin=-2, vmax=2,s=40)

ax.set_aspect('equal')
ax.set_ylim(-200,0)
ax.set_xlim(200,1000)
def dynamics_pressure3(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(phase),len(phase[0])) 
    return x,z,dypre
x,z,dypre = dynamics_pressure3(frame)
cb_plot=ax3.pcolormesh(ele_x,ele_z,dypre/1e6,cmap=cc,vmin=-200, vmax=200)#,vmin=-2, vmax=2,s=40)
ax3.set_aspect('equal')
ax3.set_ylim(-200,0)
ax3.set_xlim(200,1000)
