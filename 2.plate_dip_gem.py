#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs
from creat_database import oceanic_slab,nodes_to_elements
model = str(sys.argv[1])
#model = 'k0211'
path = '/home/jiching/geoflac/'+model+'/'
path = '/scratch2/jiching/03model/'+model+'/'
#path = '/scratch2/jiching/'+model+'/'
#path = 'F:/model/'+model+'/'
savepath = '/home/jiching/geoflac/data'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
fig, (ax)= plt.subplots(1,1,figsize=(17,12))

plot_dip = 0

phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
angle = np.zeros(end)

rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))
stslab = 0;xmean=0;ictime=20;width=300
for i in range(1,end):
    crust_x,crust_z = oceanic_slab(i)
    x, z = fl.read_mesh(i)
    ele_x, ele_z = nodes_to_elements(x,z)
    phase = fl.read_phase(i)
    trench_ind = np.argmin(z[:,0]) 
    x_trench = ele_x[:,0][np.argmin(ele_z[:,0])]
    within_plot = (ele_x[:,0]>x_trench-width)* (crust_z < 0)
    ax.plot(crust_x[within_plot]-(x_trench),crust_z[within_plot],color=newcolors[i],zorder=1)
    if not True in (crust_z < -80):
        print(i)
        continue
    if i >=end-ictime:
        stslab += crust_z
        xmean += (crust_x-x_trench)
    if plot_dip:
        ind_within_80km = (crust_z >= -80) * (crust_z < -5)
        crust_xmin = np.amin(crust_x[ind_within_80km])
        crust_xmax = np.amax(crust_x[ind_within_80km])
        crust_zmin = np.amin(crust_z[ind_within_80km])
        crust_zmax = np.amax(crust_z[ind_within_80km])
        dx = crust_xmax - crust_xmin
        dz = crust_zmax - crust_zmin
        angle[i] = math.degrees(math.atan(dz/dx))
ax.plot((xmean[within_plot]/ictime),(stslab[within_plot]/ictime),c='green',lw=3)
ax.set_aspect('equal')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Depth (km)')
ax.set_xlim(-10,450)
# ax.set_ylim(-200,0)
ax.set_title('Geometry of Subducted Slab')
fs.save_2txt(str(model)+'_stack_slab',savepath,xmean[within_plot]/ictime,stslab[within_plot]/ictime)
fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_gem.jpg')
if plot_dip:
    fig2, (ax2)= plt.subplots(1,1,figsize=(17,12))
    ax2.plot(fl.time[angle>0],angle[angle>0],c='blue',lw=2)
    ax2.set_xlim(0,fl.time[-1])
    # ax2.set_xticks(np.linspace(0,30,6))
    ax2.set_title('Angle Variation')
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Angel ($^\circ$)')
    fig2.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip.jpg')
