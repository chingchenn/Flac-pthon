#!/usr/bin/env python

import time
import flac
import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import function_for_flac as f2
import matplotlib.pyplot as plt
fig, (ax)= plt.subplots(1,1,figsize=(25,10))
start = time.time()
model = str(sys.argv[1]) 
path = '/home/jiching/geoflac/'+model+'/'
#path = '/scratch2/jiching/sem02model/'+model+'/'
#path = '/scratch2/jiching/22winter/'+model+'/'
#model='s1518'
#path = '/Volumes/My Book/model/v1'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1

cm = plt.cm.get_cmap('RdYlBu')
i=99
x,z=fl.read_mesh(i)
phase = fl.read_phase(i)
ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
pre=fl.read_pres(i)
onepre=-pre.flatten()
a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))*100
cb_plot=ax.scatter(ele_x,ele_z,c=dypre,cmap=cm,vmin=-200, vmax=200)
ax_cbin = fig.add_axes([0.63, 0.28, 0.23, 0.03])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
tick_font_size = 10
cb.ax.tick_params(labelsize=tick_font_size)
ax_cbin.set_title('MPa',fontsize=20)
ax.set_aspect('equal')
ax.set_xlim(0,max(ele_x[:,0]))
ax.set_ylim(min(ele_z[0,:]),0)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_flame='+str(i)+'_dynamic_pressure.png')
#fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')
