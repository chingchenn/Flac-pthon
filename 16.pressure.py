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
    pre = fl.read_pres(frame) *1e8
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = nodes_to_elements(x, z)
    phase = fl.read_phase(frame)
    new_pre = np.zeros((len(x)-1, len(x[0])-1))
    for xx in range(len(ele_x)):
        den = fl.read_density(frame)
        if phase[xx,0]==2:
            ref_den = den[-4,:]
            new_pre[xx,:] = pre[xx,:]-ref_den*g*ele_z[xx,:]*1e3
        else:
            ref_den = den[6,:]
            new_pre[xx,:] = pre[xx,:]-ref_den*g*ele_z[xx,:]*1e3
    return x,z,new_pre
x,z,dypre = dynamics_pressure(frame)
fig2, (ax2)= plt.subplots(1,1,figsize=(25,10))
cb_plot=ax2.pcolormesh(x,z,dypre/1e6,cmap=cc,vmin=-200, vmax=200)
ax2.set_xlim(200,1000)
ax2.set_ylim(-200,0)
ax_cbin = fig2.add_axes([0.43, 0.28, 0.23, 0.03])
cb2 = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
tick_font_size = 20
cb2.ax.tick_params(labelsize=tick_font_size)
ax_cbin.set_title('MPa',fontsize=20)
ax2.set_aspect('equal')