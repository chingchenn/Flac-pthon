#!/usr/bin/env python
import math
import time
import flac
import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import function_for_flac as f2
import matplotlib.pyplot as plt
fig, (ax)= plt.subplots(1,1,figsize=(15,15))
start = time.time()
model = str(sys.argv[1]) 
path = '/home/jiching/geoflac/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
melt = np.zeros(end)
magma = np.zeros(end)
kkmelt = np.zeros(end)
kkchamber=np.zeros(end)
rrr=np.zeros(end)
for i in range(110,111):
    mm=fl.read_pres(i)
    x,z=fl.read_mesh(i)
    ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
    ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
    melt[i]=np.max(mm)
    ax.scatter(ele_z[123,:],-mm[123,:],color='tomato')
    pre=fl.read_pres(i)-mm[123,:]
    
#ax.set_title('melt * area',fontsize=20)
#ax.set_ylim(0,0.8)
#ax.set_xlim(0,40)
phase=fl.read_phase(i)
plt.scatter(x,z,phase)

#fig2,(ax,ax2)=plt.subplots(1,2,figsize=(25,12))
#cb_plot=ax.scatter(melt,magma,c=fl.time,cmap='rainbow')
#ax_cbin = fig2.add_axes([0.13, 0.78, 0.23, 0.03])
#cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
#ax_cbin.set_title('Myr (km)')
#rrr1=f2.moving_window_smooth(rrr,12)
#ax2.plot(fl.time,rrr1,color='k',lw=3)
#ax2.plot(fl.time,rrr,color='gray',linestyle=':')
#end = time.time()
#print(end - start)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_magma_parameter_time_series.png')
#fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')
