#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:17:25 2022

@author: ji-chingchen
"""
import sys,os
import math
import numpy as np
import flac
import flacmarker2vtk
import function_for_flac as fd
import matplotlib.pyplot as plt

path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
#path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'

model = 'ch1404'
os.chdir(path+model)
fl = flac.Flac()
time=fl.time
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
bwith=3

def trench(start_vts=1,model_steps=end):
    trench_x=np.zeros(end)
    trench_z=np.zeros(end)
    trench_index=np.zeros(end)
    for i in range(start_vts,model_steps):
        x,z = fl.read_mesh(i)
        sx,sz=fd.get_topo(x,z)
        arc_ind,trench_ind=fd.find_trench_index(z)
        trench_index[i]=trench_ind
        trench_x[i]=sx[trench_ind]
        trench_z[i]=sz[trench_ind]
    return trench_index,trench_x,trench_z

start_vts=1
def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt
trench_index,trench_x,trench_z=trench()
def oceanic_slab(frame):
    phase_oceanic = 3
    phase_ecolgite = 13
    phase_oceanic_1 = 17
    phase_ecolgite_1 = 18
    x, z = fl.read_mesh(i)
    ele_x, ele_z = nodes_to_elements(x,z)
    phase = fl.read_phase(i)
    trench_ind = np.argmin(z[:,0]) 
    crust_x = np.zeros(nex)
    crust_z = np.zeros(nex)
    for j in range(trench_ind,nex):
        ind_oceanic = (phase[j,:] == phase_oceanic) + (phase[j,:] == phase_ecolgite)+(phase[j,:] == phase_oceanic_1) + (phase[j,:] == phase_ecolgite_1)
        if True in ind_oceanic:
            kk = ele_z[j,ind_oceanic]
            xx = ele_x[j,ind_oceanic]
            if len(kk[kk<-15])==0:
                continue    
            crust_x[j] = np.max(xx[kk<-15])
            crust_z[j] = np.max(kk[kk<-15])
    return crust_x,crust_z
#####===============================low viscosity area===============================
# fig,ax=plt.subplots(2,1,figsize=(10,6))
# # def read_wedgevis(trench_indexend,depth1=100, depth2=120):
# color=['r','g','b','k']
# # for www,depth1 in enumerate([80,90,100,110]):
# depth1 = 80
# depth2 = 130
# viswedge=np.zeros(end)
# areawedge=np.zeros(end)

# for i in range(1,end+1):
#     x, z = fl.read_mesh(i)
#     vis=fl.read_visc(i)
#     area=fl.read_area(i)
#     ele_x, ele_z = nodes_to_elements(x,z)
#     crust_x,crust_z = oceanic_slab(i)
#     temp = fl.read_temperature(i)
#     fig2,aa=plt.subplots(1,1,figsize=(10,6))
#     aa.scatter(ele_x,ele_z,c=vis)

#     wedge_area = np.zeros(nex-int(trench_index[i-1]))
#     for ii in range(int(trench_index[i-1]),nex):
#     # for ii in range(135,185):
#         if crust_z[ii]<-depth2:
#             break
#         up= (ele_z[ii,:]> crust_z[ii])*(ele_z[ii,:]<-depth1)*(vis[ii,:]<22)
#         if True in up:
#             wedge_area[ii-int(trench_index[i-1])]=np.mean(area[ii,up]/1e6)
#             areawedge[i-1]+=sum(area[ii,up]/1e6)
#             viswedge[i-1]+=np.mean(vis[ii,up])
#             aa.scatter(ele_x[ii,up],ele_z[ii,up],c='w',s=300)
#     if len(wedge_area[wedge_area>0])==0:
#         continue
#     areawedge[i-1] = areawedge[i-1]
#     viswedge[i-1] = viswedge[i-1]/len(wedge_area[wedge_area>0])
#     cx=aa.contour(ele_x,ele_z,vis,cmap = 'rainbow_r',levels =[23,24]) 
#     aa.contour(x,z,temp,cmap = 'magma',levels =[700])
#     aa.set_aspect('equal')
#     aa.set_xlim(400,850)
#     aa.set_ylim(-200,0)

# ax[0].scatter(fl.time[areawedge>0],areawedge[areawedge>0])
# ax[1].scatter(fl.time[areawedge>0],viswedge[areawedge>0])

# for qq in range(len(ax)):
#     ax[qq].set_xlim(0,30)
#     # ax[qq].legend()
    
    # return wedgevis
#####=====================low viscosity area (x direction)=====================
fig,ax=plt.subplots(3,1,figsize=(10,6))
# def read_wedgevis(trench_indexend,depth1=100, depth2=120):
color=['r','g','b','k']
# for www,depth1 in enumerate([80,90,100,110]):
distance_from_trench = 550
depth1 = 80
depth2 = 150
model_list=['ch1404','ch1406']
for yy,model in enumerate(model_list):
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    nex = fl.nx - 1
    nez = fl.nz - 1
    viswedge=np.zeros(end)
    areawedge=np.zeros(end)
    temwedge=np.zeros(end)
    
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        vis=fl.read_visc(i)
        area=fl.read_area(i)
        ele_x, ele_z = nodes_to_elements(x,z)
        crust_x,crust_z = oceanic_slab(i)
        temp = fl.read_temperature(i)
        ele_tem = temp_elements(temp)
        # fig2,aa=plt.subplots(1,1,figsize=(10,6))
        # aa.scatter(ele_x,ele_z,c=vis)
    
        wedge_area = np.zeros(nex-int(trench_index[i-1]))
        for ii in range(int(trench_index[i-1]),nex):
        # for ii in range(135,185):
            crust_zdonw = crust_z[ii]
            if crust_z[ii]==0 and ii>120:
                crust_zdonw = min(crust_z)
            up= (ele_z[ii,:]> -depth2)*(ele_x[ii,:]<(trench_x[i]+distance_from_trench))*(ele_z[ii,:]< -depth1)*(ele_z[ii,:]>crust_zdonw)
            if True in up:
                wedge_area[ii-int(trench_index[i-1])]=np.mean(area[ii,up]/1e6)
                areawedge[i-1]+=sum(area[ii,up]/1e6)
                viswedge[i-1]+=np.mean(vis[ii,up])
                temwedge[i-1]+=np.mean(ele_tem[ii,up])
                # aa.scatter(ele_x[ii,up],ele_z[ii,up],c='w',s=300)
        if len(wedge_area[wedge_area>0])==0:
            continue
        areawedge[i-1] = areawedge[i-1]
        viswedge[i-1] = viswedge[i-1]/len(wedge_area[wedge_area>0])
        temwedge[i-1] = temwedge[i-1]/len(wedge_area[wedge_area>0])
        # cx=aa.contour(ele_x,ele_z,vis,cmap = 'rainbow_r',levels =[23,24]) 
        # aa.contour(x,z,temp,cmap = 'magma',levels =[700])
        # aa.set_aspect('equal')
        # aa.set_xlim(200,900)
        # aa.set_ylim(-200,0)
        # aa.set_title(model+'_vis_'+str(i),fontsize=30)
        # fig2.savefig(figpath+model+'_wedge_snapshot.png')
    
    ax[0].scatter(fl.time[areawedge>0],areawedge[areawedge>0],c=color[yy],label=model)
    ax[1].scatter(fl.time[areawedge>0],viswedge[areawedge>0],c=color[yy],label=model)
    ax[2].scatter(fl.time[areawedge>0],temwedge[areawedge>0],c=color[yy],label=model)
    ax[0].set_title(model+'_vis',fontsize=30)
    ax[0].set_ylabel('area (km)',fontsize=16)
    ax[1].set_ylabel('viscosity (Pa s)',fontsize=16)
    ax[2].set_ylabel('Temperature',fontsize=16)
    ax[-1].set_xlabel('Time (Myr)',fontsize=16)
    for qq in range(len(ax)):
        ax[qq].set_xlim(0,30)
        ax[qq].tick_params(axis='x', labelsize=16)
        ax[qq].tick_params(axis='y', labelsize=16)
        ax[qq].grid()
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].legend()
fig.savefig(figpath+model_list[0]+'_'+model_list[-1]+'_wedge_compare.png')
