#!/usr/bin/env python
from locale import format_string
import math
from socketserver import ForkingTCPServer
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as fd
from Main_creat_database import oceanic_slab,nodes_to_elements
model = str(sys.argv[1])
#model = 'Ref'
path = '/home/jiching/geoflac/'+model+'/'
#path = '/scratch2/jiching/03model/'+model+'/'
#path = '/scratch2/jiching/'+model+'/'
#path = 'F:/model/'+model+'/'
#path = 'D:/model/'+model+'/'
savepath = '/home/jiching/geoflac/data'
os.chdir(path)
fl = flac.Flac();end = fl.nrec

plot_dip = 0
i = end

crust_x,crust_z = oceanic_slab(i)
x, z = fl.read_mesh(i)
ele_x, ele_z = nodes_to_elements(x,z)
phase = fl.read_phase(i)
x_trench = ele_x[:,0][np.argmin(ele_z[:,0])]
crust_xx = crust_x[crust_x>0]
crust_z = crust_z[crust_x>0]
crust_x = crust_xx
###----------------------------------dip with dix--------------------------------------
angle_distance = np.zeros(800)
ddd=np.linspace(1,800,800)
for distance_from_trench in ddd:
    ind_within = ((crust_x-x_trench) <= distance_from_trench)
    crust_xmin = np.amin(crust_x[ind_within])
    crust_xmax = np.amax(crust_x[ind_within])
    crust_zmin = np.amin(crust_z[ind_within])
    crust_zmax = np.amax(crust_z[ind_within])
    dx = crust_xmax - crust_xmin
    dz = crust_zmax - crust_zmin
    if dz == 0:
        continue
    angle_distance[int(distance_from_trench-1)] = math.degrees(math.atan(dz/dx))
fig, (ax)= plt.subplots(1,1,figsize=(14,10))
aaa = fd.moving_window_smooth(angle_distance,10)
ax.plot(ddd,aaa,lw=4,color='k')
ax.set_xlabel('Distance from the trench (km)',fontsize=20)
ax.set_ylabel('Slab dip',fontsize=20)
ax.set_xlim(0,max(crust_x-x_trench))
ax.grid()
bwith = 3
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
#ax.set_title('Geometry of Subducted Slab')
fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip_with_dix.png')
###----------------------------------dip with diz--------------------------------------

angle_depth = np.ones(300)*1
ddd=np.linspace(1,300,300)
for Depth in ddd:
    ind_within = (crust_z >= -Depth) 
    if not True in ind_within:
            continue
    crust_xmin = np.amin(crust_x[ind_within])
    crust_xmax = np.amax(crust_x[ind_within])
    crust_zmin = np.amin(crust_z[ind_within])
    crust_zmax = np.amax(crust_z[ind_within])
    dx = crust_xmax - crust_xmin
    dz = crust_zmax - crust_zmin
    if dz == 0:
        continue
    angle_depth[int(Depth-1)] = math.degrees(math.atan(dz/dx))
fig2, (ax2)= plt.subplots(1,1,figsize=(14,10))
aaa = fd.moving_window_smooth(angle_depth,10)
ax2.plot(ddd,aaa,lw=4,color='k')
ax2.set_xlabel('Depth (km)',fontsize=20)
ax2.set_ylabel('Slab dip',fontsize=20)
ax2.set_xlim(0,max(crust_x-x_trench))
ax2.grid()
bwith = 3
ax2.spines['bottom'].set_linewidth(bwith)
ax2.spines['top'].set_linewidth(bwith)
ax2.spines['right'].set_linewidth(bwith)
ax2.spines['left'].set_linewidth(bwith)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
#ax.set_title('Geometry of Subducted Slab')
fig2.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip_with_diz.png')
