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


plt.rcParams["font.family"] = "Helvetica"
#---------------------------------- DO WHAT -----------------------------------
### pdf or png
png             = 0
pdf             = 0
#---------------------------------- SETTING -----------------------------------
path = '/haome/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/22summer/'
path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
#path = 'D:/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
figpath='/Users/chingchen/Desktop/FLAC_Works/Eclogite_flat_slab/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Ref_Nazca/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/Thesis_figure/Discussion/'
# figpath='/Users/chingchen/OneDrive - 國立台灣大學/青年論壇/'
# figpath='/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/AGU/POSTER/Poster_figure/'

xmin,xmax = 250,1000
zmin,zmax = -200,10
model = 'Nazca_aa06'
frame1 = 30
frame2 = 60
frame3 = 140
frame4 = 170
limit = 3.2e-4
bwith = 4
fontsize=35
labelsize=30
skip = (slice(None, None, 8))
depth1=150
scale=0.03

os.chdir(path+model)
fl = flac.Flac()
time=fl.time


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
def oceanic_slab2(frame,x,z,phase,trench_index):
    phase_oceanic = 3
    phase_ecolgite = 13
    phase_oceanic_1 = 17
    phase_ecolgite_1 = 18
    ele_x, ele_z = flac.elem_coord(x, z)
    trench_ind = int(trench_index[frame-1])
    crust_x = np.zeros(len(ele_x))
    crust_z = np.zeros(len(ele_x))
    for j in range(trench_ind,len(ele_x)):
        ind_oceanic = (phase[j,:] == phase_oceanic) + (phase[j,:] == phase_ecolgite)+(phase[j,:] == phase_oceanic_1) + (phase[j,:] == phase_ecolgite_1)
        if True in ind_oceanic:
            kk = ele_z[j,ind_oceanic]
            xx = ele_x[j,ind_oceanic]
            if len(kk[kk<-15])==0:
                continue
            crust_x[j] = np.max(xx[kk<-15])
            crust_z[j] = np.max(kk[kk<-15])       
    return crust_x,crust_z
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
   
ax5 = ax[0] 
ax = ax[1]
frame = frame1
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
magma = fl.read_fmagma(frame)
#--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cc = plt.cm.get_cmap('jet')
cbvis=ax.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
cax = plt.axes([0.548, 0.06, 0.354, 0.02])
cc1=fig.colorbar(cbvis, ax=ax,cax=cax,orientation='horizontal',ticks=[20,21,22,23,24])
cc1.set_label(label='log$_{10}$ (viscosity) (Pa $\cdot$ s)', size=35)
cc1.ax.tick_params(labelsize=30)
cc1.ax.yaxis.set_label_position('left')
ax.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
ck = plt.cm.get_cmap('RdBu_r')
cbpre=ax5.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
cax = plt.axes([0.123, 0.06, 0.354, 0.02])
cc3=fig.colorbar(cbpre, ax=ax5,cax=cax,orientation='horizontal',ticks=[-200,-100,0,100,200])
cc3.set_label(label='dynamics pressure (MPa)', size=35)
cc3.ax.tick_params(labelsize=30)
cc3.ax.yaxis.set_label_position('left')
ax5.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax5.contour(ele_x,-ele_z,magma,colors=['#6B8E23'],levels =[limit], linewidths=6)

# xvel,zvel = fl.read_vel(frame)
# phase = fl.read_phase(frame)
# time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
# crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
# nex = len(x)
# for ii in range(int(trench_index[frame-1]),nex-1,8):
#     up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>600)
#     right=(x[ii,:]>500)*(temp[ii,:]>800)*(x[ii,:]<900)
#     if True in up :
#         ax5.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)
#     if True in right:
#         ax5.quiver(x[ii,right][skip],-z[ii,right][skip],xvel[ii,right][skip],zvel[ii,right][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)

ax5.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'k')
ax5.text(170,-5,'(a)',fontsize=fontsize+10, color = 'k')
ax.text(170,-5,'(e)',fontsize=fontsize+10, color = 'k')
#----------------------------- FIG2 -----------------------------
ax6 = ax2[0] 
ax2 = ax2[1]
frame = frame2
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
magma = fl.read_fmagma(frame)
##--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax2.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax2.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
cbpre=ax6.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
ax6.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax6.contour(ele_x,-ele_z,magma,colors=['#6B8E23'],levels =[limit], linewidths=6)

# xvel,zvel = fl.read_vel(frame)
# phase = fl.read_phase(frame)
# time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
# crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
# nex = len(x)
# for ii in range(int(trench_index[frame-1]),nex-1,8):
#     up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>600)
#     right=(x[ii,:]>600)*(temp[ii,:]>800)*(x[ii,:]<900)
#     if True in up :
#         ax6.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)
#     if True in right:
#         ax6.quiver(x[ii,right][skip],-z[ii,right][skip],xvel[ii,right][skip],zvel[ii,right][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)


ax6.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'k')
ax6.text(170,-5,'(b)',fontsize=fontsize+10, color = 'k')
ax2.text(170,-5,'(f)',fontsize=fontsize+10, color = 'k')
#----------------------------- FIG3 -----------------------------
ax7 = ax3[0] 
ax3 = ax3[1]
frame = frame3
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
magma = fl.read_fmagma(frame)
##--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax3.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax3.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)

#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
cbpre=ax7.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
ax7.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax7.contour(ele_x,-ele_z,magma,colors=['#6B8E23'],levels =[limit], linewidths=6)

# xvel,zvel = fl.read_vel(frame)
# phase = fl.read_phase(frame)
# time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
# crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
# nex = len(x)
# for ii in range(int(trench_index[frame-1]),nex-1,8):
#     up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>600)
#     if True in up :
#         ax7.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)


ax7.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'k')
ax7.text(170,-5,'(c)',fontsize=fontsize+10, color = 'k')
ax3.text(170,-5,'(g)',fontsize=fontsize+10, color = 'k')
#----------------------------- FIG4 -----------------------------
ax8 = ax4[0] 
ax4 = ax4[1]
frame = frame4
xt,zt = fl.read_mesh(frame)
ele_x,ele_z = flac.elem_coord(xt, zt)
temp = fl.read_temperature(frame)
magma = fl.read_fmagma(frame)
##--------------------- viscosity plotting -------------------------
x,z,ele_x,ele_z,vis,ztop = get_vis(frame)
cbvis=ax4.pcolormesh(ele_x,-ele_z,vis,cmap=cc,vmin=20, vmax=24,shading='gouraud')
ax4.contour(xt,-zt,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=2)
#--------------------- pressure plotting -------------------------
x,z,new_pre = dynamics_pressure(frame)
new_pre = new_pre-np.median(new_pre)
cbpre=ax8.pcolormesh(ele_x,-ele_z,new_pre/1e6,cmap=ck,vmin=-200, vmax=200,shading='gouraud')
ax8.contour(x,-z,temp,colors='#696969',levels =[200,400,600,800,1000,1200],linewidths=3)
ax8.contour(ele_x,-ele_z,magma,colors=['#6B8E23'],levels =[limit], linewidths=6)

# xvel,zvel = fl.read_vel(frame)
# phase = fl.read_phase(frame)
# time, trench_index,trench_x,trench_z = np.loadtxt(savepath+'trench_for_'+model+'.txt').T
# crust_x,crust_z = oceanic_slab2(frame, x, z, phase, trench_index)
# nex = len(x)
# for ii in range(int(trench_index[frame-1]),nex-1,8):
#     up= (z[ii,:]> crust_z[ii])*(z[ii,:]>-depth1)*(z[ii,:]<-50)*(temp[ii,:]>600)
#     if True in up :
#         vector=ax8.quiver(x[ii,up][skip],-z[ii,up][skip],xvel[ii,up][skip],zvel[ii,up][skip],
#                  scale_units='xy', scale=scale,headlength=5,headwidth=3)

ax8.text(270,170,str(np.round(fl.time[frame-1],1))+' Myr',fontsize=fontsize, color = 'k')
ax8.text(170,-5,'(d)',fontsize=fontsize+10, color = 'k')
ax4.text(170,-5,'(h)',fontsize=fontsize+10, color = 'k')
#ax8.quiverkey(vector,1070,277,3,"3 cm/y",coordinates='data',color='k',fontproperties={'size': labelsize})
# ---------------------- plot setting --------------------------
xmajor_ticks=np.array([300,400,500,600,700,800,900,1000])
ymajor_ticks = np.linspace(200,0,num=5)
for kk,aa in enumerate([ax,ax2,ax3,ax4,ax5,ax6,ax7,ax8]):
    if kk>=4:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in',pad=10)
    else:
        aa.tick_params(labelsize=labelsize,width=3,length=15,right=True, top=True,direction='in',pad=10,color='w')
    aa.set_aspect('equal')
    aa.set_yticks(ymajor_ticks)
    aa.set_ylabel('depth (km)',fontsize=fontsize)
    aa.set_xticks(xmajor_ticks)
    aa.set_xlim(xmin,xmax)
    aa.set_ylim(-zmin,-zmax)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
ax4.set_xlabel('distance (km)',fontsize=fontsize)
ax8.set_xlabel('distance (km)',fontsize=fontsize)
# ax.set_title(model,fontsize=30)
if png:
    fig.savefig(figpath+'fig4a.png')
if pdf:
    fig.savefig(figpath+'fig4av6.pdf')