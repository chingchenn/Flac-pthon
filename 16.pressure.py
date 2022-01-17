#!/usr/bin/env python

import time
import flac
import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import function_for_flac as f2
import matplotlib.pyplot as plt
fig, (ax)= plt.subplots(1,1,figsize=(25,10))
start = time.time()
model = str(sys.argv[1]) 
path = '/home/jiching/geoflac/'+model+'/'
#model='s1518'
#path = '/Volumes/My Book/model/v1'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
melt = np.zeros(end)
magma = np.zeros(end)
kkmelt = np.zeros(end)
kkchamber=np.zeros(end)
rrr=np.zeros(end)
cm = plt.cm.get_cmap('RdYlBu')
i=100
x,z=fl.read_mesh(i)
phase = fl.read_phase(i)
ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
phase=fl.read_phase(i)
pre=fl.read_pres(i)
onepre=-pre.flatten()
a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))
cb_plot=ax.scatter(ele_x,ele_z,c=dypre,vmin=-8, vmax=8,cmap=cm)
ax_cbin = fig.add_axes([0.63, 0.28, 0.23, 0.03])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')

fig2, (ax2)= plt.subplots(1,1,figsize=(25,10))
cb_plot=ax2.scatter(ele_x,ele_z,c=-pre,vmin=-8, vmax=8,cmap=cm)
ax_cbin = fig2.add_axes([0.63, 0.28, 0.23, 0.03])
cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ax.set_aspect('equal')
ax2.set_aspect('equal')
ax_cbin.set_title('kbar')

fig.savefig('/home/jiching/geoflac/figure/'+model+'_flame='+str(i)+'_dynamic_pressure.png')
#fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')
