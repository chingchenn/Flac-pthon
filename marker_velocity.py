#!/usr/bin/env python
import flac
import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_for_flac as f2

model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model
#sys.path.append('/home/jiching/geoflac/util/')
os.chdir(path)
fl = flac.Flac();end = fl.nrec
plt.rcParams['figure.figsize'] =10,8
number_of_marker = 50
xmax,xmin=760,450
zmax,zmin=-20,-50
dis = np.zeros(end)
t=np.linspace(1,end,end)
fig, (ax)= plt.subplots(1,1,figsize=(10,8))
for i in range(51,end):
    x, z, age, ph, idd, a1, a2, ntriag = fl.read_markers(i+1)
    xp,zp,agep,phase_p ,ID_p, a1, a2, ntriag = fl.read_markers(i)
    ind_mat = idd[(x>=xmin)*(x<=xmax)*(z<zmax)*(z>zmin)]
    ind_p = ID_p[(xp>=xmin)*(xp<=xmax)*(zp<zmax)*(zp>zmin)]
    marker = []
    for ind in ind_p:
        if ind in ind_mat:
            if len(marker)<=(number_of_marker-1):
                marker.append(ind)

    rrr1 = np.zeros((number_of_marker,2))
    rrr2 = np.zeros((number_of_marker,2))
    for number in range(len(marker)):
        md1=ID_p==marker[number]
        md2=idd==marker[number]
        rrr1[number]=xp[md1],zp[md1]
        rrr2[number]=x[md2],z[md2]

    x1=np.average(rrr1[:,0])
    z1=np.average(rrr1[:,1])
    x2=np.average(rrr2[:,0])
    z2=np.average(rrr2[:,1])
    dis[i]=np.sqrt((x2-x1)**2+(z2-z1)**2)
v=dis*10**6/2e5
vvv=np.array(f2.moving_window_smooth(v,8))
ax.plot(t[dis>0]*0.2,v[v>0],c='gray',linestyle=':')
ax.plot(t[dis>0]*0.2,vvv[dis>0],c='k')
print(np.average(v[v>0]),np.average(vvv[vvv>0]))
ax.set_title(str(np.average(v[v>0]))+'   '+str(np.average(vvv[vvv>0])))
ax.set_xlabel('Time (Myr)',fontsize=20)
ax.set_ylabel('velocity (mm/y)',fontsize=22)
fig.savefig('/home/jiching/geoflac/figure/'+model+'velocity.png')
