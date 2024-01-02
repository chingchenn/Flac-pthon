#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 22:51:52 2023

@author: chingchen
"""


import os
import flac
import numpy as np
import matplotlib
import function_for_flac as fd
import matplotlib.pyplot as plt
from P1fig5b import oceanic_slab2
#import flac_interpolate as fi


plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- DO WHAT -----------------------------------
### pdf or png
png             = 0
pdf             = 1

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath = '/Users/chingchen/Desktop/Flac_Works/Eclogite_flat_slab/'


model1 = 'Nazca_aa06'
model2 = 'Nazca_aa06'
model2 = 'Nazca_a0706'
frame1 = 55

limit = 1e-4
skip = (slice(None, None, 3))
depth1=150
scale=0.15
headlength=5
headwidth=3

bwith = 4
fontsize=30
labelsize=30
shading = 'gouraud'
#shading = 'nearest'
#------------------------------------------------------------------------------
def plot_snapshot(frame):
    x,z = fl.read_mesh(frame)
    xtop,ztop = fd.get_topo(x,z)
    phase = fl.read_phase(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    temp = fl.read_temperature(frame)
    return x, z, ele_x, ele_z, phase, temp, ztop
def get_vis(frame):
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    vis = fl.read_visc(frame)
    xtop,ztop = fd.get_topo(x,z)
    return x,z,ele_x,ele_z,vis,ztop
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
      "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
      "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
      "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
colors2=[
 '#C98F49', '#92C0DF', '#2553C7', '#FFFFFF', '#6495ED',
 '#2E8B57', '#524B52', '#9A32CD', '#6B8E23','#D4DBF4',
 '#D8BFD8','#999999','#F2C85B','#92C0DF','#999999',
 '#4CC552','#999999','#999999','#999999','#999999']
phase8= matplotlib.colors.ListedColormap(colors2)

def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre
#------------------------------------------------------------------------------
#fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,22))
fig, (ax,ax2)= plt.subplots(2,2,figsize=(23,12))
xmin,xmax = 250,1000
zmin,zmax = -300,10
#----------------------------- FIG1 -----------------------------
model = model1
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
limit = 1e-3
#----------------------------- FIG 1-1 -----------------------------    
ax5 = ax2[0]
ax6 = ax2[1] 
ax2 = ax[1] 
ax = ax[0]
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
ck = plt.cm.get_cmap('RdBu_r')
cbpre=ax.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading=shading)
cax3 = plt.axes([0.123, 0.04, 0.354, 0.02])
cc3=fig.colorbar(cbpre, ax=ax,cax=cax3,orientation='horizontal',ticks=[-200,-100,0,100,200])
cc3.set_label(label='dynamics pressure (MPa)', size=35)
cc3.ax.tick_params(labelsize=30)
cc3.ax.yaxis.set_label_position('left')
ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax.text(270,270,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
#----------------------------- FIG 1-2 -----------------------------    
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax2.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading=shading)
cax1 = plt.axes([0.548, 0.04, 0.354, 0.02])
cc1=fig.colorbar(cbvis, ax=ax,cax=cax1,orientation='horizontal',ticks=[20,21,22,23,24])
cc1.set_label(label='viscosity (Pa$\cdot$s)', size=35)
cc1.ax.tick_params(labelsize=30)
cc1.ax.yaxis.set_label_position('left')
ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG 2 -----------------------------
model = model2
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
#limit = 1e-4
#----------------------------- FIG 2-1-----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#magma = fl.read_fmagma(frame)
#magma_limit = (magma>limit)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
cbpre=ax5.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading=shading)
ax5.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG 2-2-----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax6.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading=shading)
ax6.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax5.text(270,270,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize)
# ---------------------- plot setting --------------------------
xmajor_ticks=np.array([250,400,600,800,1000])

ymajor_ticks = np.linspace(300,0,num=4)
for aa in [ax,ax2,ax5,ax6]:
    if aa==ax or aa==ax5:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in',pad=13)
    else:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in',pad=13,color='w')

    aa.set_aspect('equal')
    aa.set_yticks(ymajor_ticks)
    aa.set_ylabel('depth (km)',fontsize=fontsize)
    aa.set_xticks(xmajor_ticks)
    aa.set_xlim(xmin,xmax)
    aa.set_ylim(-zmin,-zmax)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
ax5.set_xlabel('distance (km)',fontsize=fontsize)
ax6.set_xlabel('distance (km)',fontsize=fontsize)
ax.set_title('h = 5 km',fontsize=35)
ax5.set_title('h=7.5 km',fontsize=35)
############################## ZOOM IN FIGURE##############################

fig2, (ax11,ax12)= plt.subplots(2,2,figsize=(18,12))
xmin,xmax = 470,620
zmin,zmax = -150,-50
#----------------------------- FIG1 -----------------------------
model = model1
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
limit = 1e-3
nex = len(x)
#----------------------------- FIG 1-1 -----------------------------    
ax15 = ax12[0]
ax16 = ax12[1] 
ax12 = ax11[1] 
ax11 = ax11[0]
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
ck = plt.cm.get_cmap('RdBu_r')
cbpre=ax11.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading=shading)
ax11.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#-------------------- velocity vector ---------------------------
xvel,zvel = fl.read_vel(frame)
phase = fl.read_phase(frame)
time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
for ii in range(int(trench_index[frame-1]),nex-1,3):
    up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>500)
    right=(x[ii,:]>500)*(temp[ii,:]>800)*(x[ii,:]<700)
    if True in up :
        ax11.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
                 scale_units='xy', scale=scale,headlength=headlength,headwidth=headwidth)
    if True in right:
        ax11.quiver(x[ii,right][skip],-z[ii,right][skip],xvel[ii,right][skip],zvel[ii,right][skip],
                 scale_units='xy', scale=scale,headlength=headlength,headwidth=headwidth)
#----------------------------- FIG 1-2 -----------------------------    

frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax12.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading=shading)
ax12.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG 2 -----------------------------
model = model2
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
#----------------------------- FIG 2-1-----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
cbpre=ax15.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading=shading)
ax15.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#-------------------- velocity vector ---------------------------
xvel,zvel = fl.read_vel(frame)
phase = fl.read_phase(frame)
time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
for ii in range(int(trench_index[frame-1]),nex-1,5):
    up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>500)
    right=(x[ii,:]>500)*(temp[ii,:]>800)*(x[ii,:]<700)
    if True in up :
        vector=ax15.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
                  scale_units='xy', scale=scale,headlength=headlength,headwidth=headwidth)
    if True in right:
        ax15.quiver(x[ii,right][skip],-z[ii,right][skip],xvel[ii,right][skip],zvel[ii,right][skip],
                  scale_units='xy', scale=scale,headlength=headlength,headwidth=headwidth)
ax15.quiverkey(vector,540,175,5,"10 cm/y",coordinates='data',color='k',fontproperties={'size': labelsize})

#----------------------------- FIG 2-2-----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax16.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading=shading)
ax16.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
# ---------------------- plot setting --------------------------
xmajor_ticks=np.array([470,520,570,620])

ymajor_ticks = np.array([50.,100.,150.])
for kk,aa in enumerate([ax11,ax15,ax12,ax16]):
    aa.tick_params(axis = 'x',pad=15)
    aa.set_aspect('equal')
    aa.set_yticks(ymajor_ticks)
    aa.set_xticks(xmajor_ticks)
    aa.set_xlim(xmin,xmax)
    aa.set_ylim(-zmin,-zmax)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    if kk < 2:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
    else:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in',color='w')
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_color('white')
    if kk==3:
        aa.set_xlabel('distance (km)',fontsize=35)

#fig2.savefig(figpath+'fig5a_zoomin_v3.pdf')
if png:
    fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure_compare.png')
if pdf:
    #fig.savefig(figpath+'fig5a_v4.pdf')
    #fig2.savefig(figpath+'fig5a_v4zoomin.pdf')
    fig.savefig(figpath+'fig6c_v2.pdf')
    #fig2.savefig(figpath+'fig6c_v1zoomin.pdf')