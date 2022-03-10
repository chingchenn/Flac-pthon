#!/usr/bin/env python
# 2022 feb.23
import flac
import sys,os
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs
from creat_database import trench, get_topo, nodes_to_elements

model = sys.argv[1]
#model = 'h0133'
path = '/home/jiching/geoflac/'+model+'/'
#path='/scratch2/jiching/03model/'+model+'/'
#path = 'F:/model/'+model+'/'
os.chdir(path)
    
fl = flac.Flac()
end = fl.nrec

rainbow = cm.get_cmap('gray_r',end)    
newcolors = rainbow(np.linspace(0, 1, end))
cmap = plt.cm.get_cmap('gist_earth')
zmax, zmin =10, -10
xx2 = np.zeros(end) 
zz2 = np.zeros(end)
arc_x = np.zeros(end)
savepath = '/home/jiching/geoflac/data/'
trenchfile='/home/jiching/geoflac/'+'data/trench_for_'+model+'.csv'
dis, time, topo = get_topo()
if not os.path.exists(trenchfile):
    print('No file')
    trench_index,trench_x,trench_z = trench()
    name='trench_for_'+model
    fs.save_3array(name,savepath,time,trench_x,trench_z,'time','trench_x','trench_z')

df = pd.read_csv(trenchfile)
width = 600
fig, (ax) = plt.subplots(1,1,figsize=(10,12))
qqq=ax.scatter(dis,time,c=topo,cmap='gist_earth',vmax=6,vmin=-10)
cbar=fig.colorbar(qqq,ax=ax)
ax.plot(df.trench_x[df.trench_x>0],df.time[df.trench_x>0],c='k',lw=2)
ax.set_xlim(0,dis[-1][-1])
ax.set_ylim(0,fl.time[-1])
ax.set_title(str(model)+" Bathymetry Evolution",fontsize=24)
ax.set_ylabel('Time (Myr)',fontsize=20)
ax.set_xlabel('Distance (km)',fontsize=20)
cbar.set_label('Topography (km)',fontsize=20)
fig.savefig('/home/jiching/geoflac/figure/'+model+'_topo.png')

fig2, (ax2) = plt.subplots(1,1,figsize=(10,6))
topo1 = 0; aa=0
xmean = 0
ictime = 20
for i in range(1,end):
    x, z = fl.read_mesh(i)
    xt = x[:,0]
    zt = z[:,0]
    t = np.zeros(xt.shape)
    t[:] = i*0.2
    # print(xt[np.argmin(zt)])
    arc_x[i] = xt[np.argmax(zt)]
    x_mid = df.trench_x[i]
    within_plot = (xt>x_mid-width) * (xt<x_mid+width)
    if i >= end-ictime:
        topo1 += zt
        xmean += (xt-x_mid)
    ax2.plot(xt[within_plot]-x_mid,zt[within_plot],c=newcolors[i])
    

#ind_within=(arc_x<900)*(trench_x>100)

ax2.plot((xmean[within_plot]/ictime), (topo1[within_plot]/ictime),c='purple')
ax2.set_xlim(-width,width)
ax2.set_title(str(model)+"Topography")
ax2.set_ylabel("Bathymetry (m)")
ax2.set_xlabel("Distance relative to trench (km)")
fs.save_2txt(str(model)+'_stack_topography',savepath,xmean[within_plot]/ictime,topo1[within_plot]/ictime)

#ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
#cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
#cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
#ax_cbin.set_title('Bathymetry (km)')
fig2.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo_analysis.jpg')
