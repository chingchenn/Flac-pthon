#!/usr/bin/env python

import flac
import os,sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import function_for_flac as f2
import matplotlib.pyplot as plt
fig, (ax)= plt.subplots(1,1,figsize=(25,10))
#model = sys.argv[1]
model = 'b0702k'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Users/ji-chingchen/Desktop/data/'
savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1

cm = plt.cm.get_cmap('RdYlBu_r')
i=55
x,z=fl.read_mesh(i)
phase = fl.read_phase(i)
ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
pre=fl.read_pres(i)*1e8 #kb = 1e3 bar = 1e3*1e5 Pa = 1e8 N/m^2
onepre=-pre.flatten()
a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))/1e6 #MPa
cb_plot=ax.scatter(ele_x,ele_z,c=-dypre,cmap=cm)#,vmin=-2, vmax=2,s=40)
ax_cbin = fig.add_axes([0.43, 0.28, 0.23, 0.03])
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
tick_font_size = 20
cb.ax.tick_params(labelsize=tick_font_size)
ax_cbin.set_title('MPa',fontsize=20)
ax.set_aspect('equal')
ax.set_xlim(0,max(ele_x[:,0]))
ax.set_ylim(min(ele_z[0,:]),0)
ax.set_xlim(500,700)
ax.set_ylim(-200,0)
#fig.savefig('/home/jiching/geoflac/figure/'+model+'_flame='+str(i)+'_dynamic_pressure.png')
#fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')


fig2, (ax2)= plt.subplots(1,1,figsize=(25,10))
density = fl.read_density(i)
pressure = fl.read_pres(i)
static = (-density.flatten()*10*1000*ele_z.flatten()).reshape(len(phase),len(phase[0])) # N/m^2
dypre_n = (-pressure*1e8-static)/1e6 #MPa

cb_plot=ax2.scatter(ele_x,ele_z,c=-dypre_n,cmap=cm,vmin=-200, vmax=200,s=40)
ax2.set_xlim(0,max(ele_x[:,0]))
ax2.set_ylim(min(ele_z[0,:]),0)
ax2.set_xlim(500,700)
ax2.set_ylim(-200,0)
ax_cbin = fig2.add_axes([0.43, 0.28, 0.23, 0.03])
cb2 = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
tick_font_size = 20
cb2.ax.tick_params(labelsize=tick_font_size)
ax_cbin.set_title('MPa',fontsize=20)
ax2.set_aspect('equal')