#!/usr/bin/env python
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
import matplotlib
matplotlib.use('Agg')
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
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = nodes_to_elements(x,z)
    phase = fl.read_phase(frame)
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
fig3,(ax3,ax4)=plt.subplots(2,1,figsize=(20,18))
ax3.grid()
ax4.grid()
color=['#2F4F4F','#4682B4','#CD5C5C','#708090',
      '#AE6378','#282130','#7E9680','#24788F',
      '#849DAB','#EA5E51','#35838D','#4198B9',
      '#414F67','#97795D','#6B0D47','#A80359',
      '#52254F','r'] 
model_list=['ch1519','ch1522','ch1406','ch1512','ch1513','ch1510',
            'ch1520','ch1528','ch1529','ch1521','ch1530','ch1531',
            'ch1516','ch1517','ch1404','ch1523','ch1532','ch1533']
model_list=['ch1522','ch1512','ch1513','ch1510','ch1528','ch1529','ch1521',
            'ch1530','ch1531','ch1517','ch1404','ch1532','ch1519','ch1406',
            'ch1520','ch1516','ch1523','ch1533','ch1520']

depth1 = 80
depth2 = 130
for www,model in enumerate(model_list):
    os.chdir(path+model)
    fl = flac.Flac()
    time=fl.time
    end = fl.nrec
    nex = fl.nx - 1
    nez = fl.nz - 1
    viswedge=np.zeros(end)
    areawedge=np.zeros(end)
    temwedge=np.zeros(end)
    channel = np.zeros(end)
    for i in range(1,end+1):
        x, z = fl.read_mesh(i)
        vis=fl.read_visc(i)
        area=fl.read_area(i)
        ele_x, ele_z = nodes_to_elements(x,z)
        crust_x,crust_z = oceanic_slab(i)
        temp = fl.read_temperature(i)
        ele_tem = temp_elements(temp)
        wedge_area = np.zeros(nex-int(trench_index[i-1]))
        for ii in range(int(trench_index[i-1]),nex):
            if crust_z[ii]<-depth2:
                break
            up= (ele_z[ii,:]> crust_z[ii])*(ele_z[ii,:]<-depth1)*(vis[ii,:]<21)
            if True in up:
                wedge_area[ii-int(trench_index[i-1])]=np.mean(area[ii,up]/1e6)
                areawedge[i-1]+=sum(area[ii,up]/1e6)
                viswedge[i-1]+=np.mean(vis[ii,up])
                temwedge[i-1]+=np.mean(ele_tem[ii,up])
                channel[i-1]=(max(ele_z[ii,up])-min(ele_z[ii,up]))
        if len(wedge_area[wedge_area>0])==0:
            continue
        temwedge[i-1] = temwedge[i-1]/len(wedge_area[wedge_area>0])
    ccc = fd.moving_window_smooth(channel[channel>0],10)
    if www >11:
        ax3.plot(fl.time[channel>0], ccc,c = color[3],lw=3,label = model)
        ax4.plot(fl.time[channel>0], temwedge[channel>0],c = color[3],lw=3,label = model)
    else:
        ax3.plot(fl.time[channel>0], ccc,c = color[2],lw=3,label = model)
        ax4.plot(fl.time[channel>0], temwedge[channel>0],c = color[2],lw=3,label = model)
#ax3.legend(fontsize = 25)
ax3.set_ylim(0,40)
ax3.set_xlim(0,20)
#ax4.set_ylim(0,40)
ax4.set_xlim(0,20)
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.spines['bottom'].set_linewidth(bwith)
ax3.spines['top'].set_linewidth(bwith)
ax3.spines['right'].set_linewidth(bwith)
ax3.spines['left'].set_linewidth(bwith)
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax4.spines['bottom'].set_linewidth(bwith)
ax4.spines['top'].set_linewidth(bwith)
ax4.spines['right'].set_linewidth(bwith)
ax4.spines['left'].set_linewidth(bwith)
fig3.savefig('/home/jiching/geoflac/figure/'+model_list[0]+'_'+model_list[-1]+'_wedgechannel4_21.png')
