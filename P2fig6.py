#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 16:10:17 2024

@author: chingchen
"""


import flac
import os
import numpy as np
import matplotlib
import function_for_flac as fd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation/'

xmin,xmax = 250,1000
zmin,zmax = -200,10
model = 'Nazca_aa06'
model = 'Nazca_v2_06'

os.chdir(path+model)
fl = flac.Flac()
time=fl.time
bwith = 2
fontsize=33
labelsize=30
magnitude_limit = 50
#------------------------------------------------------------------------------
def compute_s1(sxx, szz, sxz):
    mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
    theta = 0.5 * np.arctan2(2*sxz,  sxx-szz)
    nx, nz = sxx.shape
    tmp = np.zeros((nx, nz, 3), dtype=sxx.dtype)
    tmp[:,:,0] = mag * np.sin(theta)
    tmp[:,:,1] = mag * np.cos(theta)
    return tmp

shift = 320
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
frame = 150
skip = (slice(None, None, 4), slice(None, None, 3))
# skip = (slice(None, None, 1), slice(None, None, 1))
fig, (ax)= plt.subplots(1,1,figsize=(17,16))
temp = fl.read_temperature(frame)
time,ele_trench,x_trench,z_trench=np.loadtxt(savepath+'trench_for_'+model+'.txt').T
chamber_limit = 1e-3
#--------------------- plotting -------------------------
### tpopgraphy
x,z = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(x, z)
xtt = x[:,0]
ztt = z[:,0]
trench_x=xtt[np.argmin(ztt)]
ele_xt = ele_x[:,0]
ele_zt = ele_z[:,0]
from sklearn.linear_model import LinearRegression
X = ele_xt
X = np.reshape(X, (len(X), 1))
y = ele_zt
modely = LinearRegression()
modely.fit(X, y)
trend = modely.predict(X)

ax.plot(ele_x[:,0]-trench_x,-(2*(ele_z[:,0]-trend)),c = 'k',lw=3,linestyle='dashed')

### stress II
sxx = fl.read_sxx(frame)*100
sxz = fl.read_sxz(frame)*100
szz = fl.read_szz(frame)*100
sII = fl.read_sII(frame)
s1,s3,s2 = compute_s1(sxx, szz, sxz).T
cbsxx = plt.cm.get_cmap('hot_r')
cbsxx=ax.pcolormesh(ele_x-trench_x,-ele_z,sII*100,cmap=cbsxx,vmin=0,vmax=1000,shading='gouraud')
magnitude = np.sqrt(s1.T[skip]**2+s3.T[skip]**2)
vector = ax.quiver(ele_x[skip][(magnitude>magnitude_limit)]-trench_x,-ele_z[skip][(magnitude>magnitude_limit)],
          s1.T[skip][(magnitude>magnitude_limit)],s3.T[skip][(magnitude>magnitude_limit)],
          pivot = 'mid',headwidth=0,color = 'green',
            scale_units='xy', scale=8,)
ax.quiverkey(vector,50,170,500,"500 MPa",coordinates='data',color='green',fontproperties={'size': labelsize-3})
cax = plt.axes([0.94, 0.340, 0.02, 0.31])
cc1=fig.colorbar(cbsxx, ax=ax,cax=cax)
cc1.ax.tick_params(labelsize=20)
cc1.set_label(label='sII (MPa)', size=25)
cc1.ax.yaxis.set_label_position('left')
ax.contour(x-trench_x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=1.5)
# ---------------------- plot setting --------------------------
ax.set_aspect('equal')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(bwith)
ax.tick_params(axis='x', labelsize=fontsize)
ax.tick_params(axis='y', labelsize=fontsize)
ymajor_ticks = np.linspace(300,0,num=7)
ax.set_yticks(ymajor_ticks)

ax.set_ylim(200,-10)
ax.set_xlim(0,700)
# ax.set_ylim(180,130)
# ax.set_xlim(550,750)

xx,zz = np.loadtxt(savepath+model+'_final_slab.txt').T
xx,zz,xt = np.loadtxt(savepath+model+'_30.0_final_slab.txt').T 
xx=xx[zz<0]
zz=zz[zz<0]
zz = fd.moving_window_smooth(zz,3)
# ax.plot(xx,-zz,color='k',lw=3)
ax.scatter(xx,-zz,color='k',s=10)
   
xxm,zzm,xtm = np.loadtxt(savepath+model+'_30_final_moho_slab.txt').T
xxm=xxm[zzm<0]
zzm=zzm[zzm<0]
zzm = fd.moving_window_smooth(zzm,3)
# ax.plot(xxm,-zzm+2,color='k',lw=3)
ax.scatter(xxm,-zzm+2,color='k',s=10)
# 
ax.set_xlabel('distance (km)',fontsize=fontsize)
ax.set_ylabel('depth (km)',fontsize=fontsize)
# fig.savefig(figpath+'figure3b_v3.pdf')