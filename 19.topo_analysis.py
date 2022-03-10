#!/usr/bin/env python
# 2022 feb.23
import flac
import sys,os
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_savedata as fs
from creat_database import trench, get_topo, nodes_to_elements

model = sys.argv[1]
#model = 'h0133'
path = '/home/jiching/geoflac/'+model+'/'
path='/scratch2/jiching/03model/'+model+'/'
#path = 'F:/model/'+model+'/'
os.chdir(path)
    
fl = flac.Flac()
end = fl.nrec
    
cmap = plt.cm.get_cmap('gist_earth')
zmax, zmin =10, -10
trench_x = np.zeros(end)
trench_t = np.zeros(end)
xx2 = np.zeros(end) 
zz2 = np.zeros(end)
arc_x = np.zeros(end)
savepath = '/home/jiching/geoflac/data/'
trenchfile='/home/jiching/geoflac/'+'data/trench_for_'+model+'.csv'
dis, time, topo = get_topo()
if not os.path.exists(trenchfile):
    trench_index,trench_x,trench_z = trench()
    name='trench_for_'+model
    fs.save_3array(name,savepath,time,trench_x,trench_z,'time','trench_x','trench_z')

df = pd.read_csv(trenchfile)

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

for i in range(end):
    x, z = fl.read_mesh(i)
    xmax = np.amax(x)
    xmin = np.amin(x)
    
    xt = x[:,0]
    zt = z[:,0]
    t = np.zeros(xt.shape)
    t[:] = i*0.2
    trench_t[i] = t[0]
    trench_x[i] = xt[np.argmin(zt)]
    # print(xt[np.argmin(zt)])
    arc_x[i] = xt[np.argmax(zt)]
    x_mid = df.trench_x[i]
    
    zz = zt[xt>xt[np.argmin(zt)]]
    z2 = zt[np.argmin(zt)+np.argmin(zz)]
    xx2[i] = xt[np.argmin(zt)+np.argmin(zz)]
    print(np.argmin(zt),np.argmin(zt)+np.argmin(zz))
    print(np.argmin(zz))
    # xx2[i] = x2

ind_within=(arc_x<900)*(trench_x>100)
#ax.plot(arc_x[ind_within],trench_t[ind_within],c='r',lw=4)
ax.plot(trench_x[ind_within],trench_t[ind_within],'k-',lw='4')
ax.plot(xx2[ind_within],trench_t[ind_within],c = 'r',lw='4')

ax.set_xlim(xmin,xmax)
ax.set_ylim(0,t[0])
ax.set_title(str(model)+" Bathymetry Evolution")
ax.set_ylabel("Time (Ma)")
ax.set_xlabel("Distance (km)")
distance=arc_x-trench_x
#ax2.plot(distance[ind_within],trench_t[ind_within],c='b',lw=4)    
#ax2.set_ylim(0,t[0])
#ax2.set_xlabel('Distance between arc and trench (km)')

ax_cbin = fig.add_axes([0.67, 0.18, 0.23, 0.03])
cb_plot = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=cmap,vmin=zmin, vmax=zmax)
cb = fig.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ax_cbin.set_title('Bathymetry (km)')
fig.savefig('/home/jiching/geoflac/figure'+'/'+str(model)+'_topo.jpg')
