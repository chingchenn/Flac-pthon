#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 19:29:35 2023

@author: chingchen
"""

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
#import flac_interpolate as fi


plt.rcParams["font.family"] = "Times New Roman"
#---------------------------------- DO WHAT -----------------------------------
### pdf or png
png             = 0
pdf             = 0

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
figpath = '/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'


xmin,xmax = 250,1000
zmin,zmax = -200,10
model1 = 'Nazca_aa16'
#model1 = 'Nazca_aa11'
model2 = 'Nazca_aa06'
#model2 = 'Nazca_aa14'
frame1 = 30
frame2 = 50
frame3 = 57
frame4 = 130

bwith = 4
fontsize=35
labelsize=30

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
fig, (ax,ax2,ax3,ax4)= plt.subplots(4,2,figsize=(34,22))
# fig, (ax,ax2,ax3)= plt.subplots(3,1,figsize=(17,17))

#----------------------------- FIG1 -----------------------------
model = model1
os.chdir(path+model)
fl = flac.Flac()
time=fl.time

#----------------------------- FIG1-1 -----------------------------    
ax5 = ax[1] 
ax = ax[0]
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
cax = plt.axes([0.945, 0.285, 0.01, 0.431])
cc3=fig.colorbar(cbvis, ax=ax,cax=cax)
cc3.set_label(label='Viscosity (Pa$\cdot$s)', size=35)
cc3.ax.tick_params(labelsize=30)
cc3.ax.yaxis.set_label_position('left')
ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'white')
#----------------------------- FIG1-2 -----------------------------    
ax6 = ax2[1] 
ax2 = ax2[0]
frame = frame2
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax2.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax2.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'white')

#----------------------------- FIG1-3 -----------------------------
ax7 = ax3[1] 
ax3 = ax3[0]
frame = frame3
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)

#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax3.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax3.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax3.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'white')

#----------------------------- FIG1-4 -----------------------------
ax8 = ax4[1] 
ax4 = ax4[0]
frame = frame4
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)

#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax4.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax4.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax4.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'white')
#----------------------------- FIG2 -----------------------------
model = model2
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
#----------------------------- FIG1-5-----------------------------
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax5.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax5.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG1-6-----------------------------
frame = frame2
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax6.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax6.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG1-7 -----------------------------
frame = frame3
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax7.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax7.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
#----------------------------- FIG1-8 -----------------------------
frame = frame4
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
#--------------------- pressure plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax8.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax8.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
# ---------------------- plot setting --------------------------
xmajor_ticks=np.array([250,300,400,500,600,700,800,900,1000])

ymajor_ticks = np.linspace(200,0,num=5)
for aa in [ax,ax2,ax3,ax4,ax5,ax6,ax7,ax8]:
    aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in')
    aa.set_aspect('equal')
    aa.set_yticks(ymajor_ticks)
    aa.set_ylabel('Depth (km)',fontsize=fontsize)
    aa.set_xticks(xmajor_ticks)
    aa.set_xlim(xmin,xmax)
    aa.set_ylim(-zmin,-zmax)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
ax4.set_xlabel('Distance (km)',fontsize=fontsize)
ax8.set_xlabel('Distance (km)',fontsize=fontsize)
# ax.set_title(model,fontsize=30)
if png:
    fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure_compare.png')
if pdf:
    fig.savefig(figpath+model+'frame_'+str(frame)+'_pressure_comparepdf')