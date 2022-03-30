#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:28:16 2021

@author: jiching
"""
import math
import flac
import os,sys
import numpy as np
from scipy import interpolate
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_for_flac as f2
## -----------------------------------------------------------
'''
Shows the oceainc crust derivative directly, smoothing, spline 
or polynimal ways. 
only for testing, too long and waste time.
One thing for sure : the Spline calculation can not used in this case.
'''
## -----------------------------------------------------------
#=========================setting=============================
model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model+'/'
#model='w1261'
#model='s1401'
path = 'D:/model/'+model+'/'
#path = '/scratch2/jiching/sem02model/'+model+'/'
#path = '/scratch/jiching/summer2021/week11/'+model+'/'
#path = '/scratch2/jiching/'+model+'/'
#path = '/Volumes/My Book/model/'+model+'/'
# path = '/Volumes/SSD500/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
#=========================Parameters=========================
phase_oceanic = 3
phase_ecolgite = 13
angle = np.zeros(end)
bet = 10
find_flat_dz1=[]
find_flat_dz2=[]
figg=0
figg2=1
fig_spline=0
mindepth=-300
#=========================main code===========================
i = 230
#for i in range(1,end):
x, z = fl.read_mesh(i)
mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
## In this code, we considered the marker phase, not the element phase
trench_ind = np.argmin(z[:,0]) 
x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]
m=[]; m2=[]
x_ocean = mx[(phase==phase_ecolgite)+(phase==phase_oceanic)]
z_ocean = mz[(phase==phase_ecolgite)+(phase==phase_oceanic)]
start = math.floor(x_trench)
final = math.floor(np.max(x_ocean))
x_grid = np.arange(start,final,bet)
ox = np.zeros(len(x_grid))
oz = np.zeros(len(x_grid))
px = start-bet
#find initial basalt depth to remove the weage basalt
kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
x_ocean = x_ocean[z_ocean<kk]
z_ocean = z_ocean[z_ocean<kk]
# interplate to the grid length "bet"
for yy,xx in enumerate(x_grid):
    oz[yy] = np.average(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
    ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
    px = xx
## ---------------------------------------------------------------
###===========================smoothing===========================
## ---------------------------------------------------------------
kkx=(f2.moving_window_smooth(ox,5))[1:-5]
kkz=(f2.moving_window_smooth(oz,5))[1:-5]
kkz=(f2.moving_window_smooth(kkz,5))[1:]
kkx=kkx[1:]
#if len(ox[oz>mindepth])<30:
#	    continue
## ------------------------ first derivative -----------------------
for kk in range(1,len(kkx)):
    cx1=kkx[kk-1];cx2=kkx[kk]
    cz1=kkz[kk-1];cz2=kkz[kk]
    if (cx2-cx1) != 0:
      m.append((cz2-cz1)/(cx2-cx1))
qq = kkx[1:]
## ------------------------ second derivative ------------------------
for ww in range(1,len(m)):
    cz1=m[ww-1];cz2=m[ww]
    cx1=qq[ww-1];cx2=qq[ww]
    if (cx2-cx1) != 0:
        m2.append((cz2-cz1)/(cx2-cx1))
    else: www = ww
qq2=qq[0:ww]
## ---------------------------- smoothing ----------------------------
mmm=f2.moving_window_smooth(m,5)
mmm2=f2.moving_window_smooth(m2,6)
### =========================== polynomial ===========================
ox = ox[oz>mindepth]
oz = oz[oz>mindepth]
z1=np.polyfit(ox,oz,4)
p4=np.poly1d(z1)
w1=p4(ox)
p3=np.polyder(p4,1)
p2=np.polyder(p4,2)
w2=p3(ox)
w3=p2(ox)
### =========================== Spline ===========================
kke=3;ss=50
new_xx=np.linspace(min(ox),max(ox),len(ox))
tck = interpolate.splrep(new_xx,oz,k=kke,s=ss)
zz0=interpolate.splev(new_xx,tck,der=0)    
zz1=interpolate.splev(new_xx,tck,der=1)    
zz2=interpolate.splev(new_xx,tck,der=2) 
yders = interpolate.spalde(ox, tck)
if fig_spline:
    fig0,(aa1,aa2,aa3)= plt.subplots(3,1,figsize=(9,12))
    aa1.plot(new_xx,zz0,c='r')
    aa1.plot(ox,oz,'k--')
    aa1.set_ylim(mindepth,0)
    aa1.grid();aa2.grid();aa3.grid()
    aa2.plot(ox,zz1)
    aa3.plot(ox,zz2)
    aa1.tick_params(axis='x', labelsize=16)
    aa1.tick_params(axis='y', labelsize=16)
    aa3.set_ylim(-0.1,0.1)
    fig0.savefig(path+str(model)+'_'+str(i)+'_slab_spline.png')
mm=-1;ff2=[]
for pp,uu in enumerate(zz2):
    if mm*uu<0:
        ff2.append(ox[pp])
    mm = uu
if len(ff2)>1 and (ff2[1]-ff2[0])>100 and ff2[0]>start:
    find_flat_dz2.append(fl.time[i])
#    if len(ff1)>1 and (ff1[1]-start)>50:
#        find_flat_dz1.append(i)
#===========================================================================
if figg:
    ### This figure show the origin point derivative directly
    fig, (bbb,aaa,ccc)= plt.subplots(3,1,figsize=(9,12)) 
    bbb.set_xlim(start,final)
    bbb.scatter(x_ocean,z_ocean,color='orange',s=20,label ='marker points') 
    bbb.scatter(ox,oz,color='cyan',s=10,label ='average marker points') 
    aaa.plot(qq,m,color='gray',zorder=1,label ='first derivative')
    aaa.plot(qq,mmm,color='k',zorder=1,label ='first derivative smoothing')
    ccc.plot([start,final],[0,0],'--',zorder=0,color='red')
    ccc.plot(qq2,m2,color='gray',label ='second derivative')
    ccc.plot(qq2,mmm2,color='k',label ='second derivative smoothing')
    bbb.set_title('frame='+str(i))
    aaa.set_xlim(start,final)
    aaa.grid();ccc.grid();bbb.grid()
    aaa.tick_params(axis='x', labelsize=16)
    aaa.tick_params(axis='y', labelsize=16)
    bbb.tick_params(axis='x', labelsize=16)
    bbb.tick_params(axis='y', labelsize=16)
    ccc.tick_params(axis='y', labelsize=16)
    ccc.set_xlim(start,final)
    ccc.set_xlabel('Distance (km)')
    fig.savefig(path+model+'frame='+str(i)+'_fig1.png')
if figg2:
    ### This figure show the origin point fitting with polynomail function
    fig2,(q1,q2,q3)= plt.subplots(3,1,figsize=(9,12))
    q1.plot(ox,w1,c='k',lw=3,label ='4st polynomial')
    q1.scatter(ox,oz,c='cyan',s=20,label ='average marker points')           
    q2.plot(ox,w2,c='k',label ='3st polynomial')
    q3.plot(ox,w3,c='k',label ='2st polynomial')
    q3.plot([start,final],[0,0],'--',zorder=0,color='red')
    q1.set_title('frame='+str(i))
    q2.set_xlim(start,final)
    q1.set_xlim(start,final) 
    q3.set_xlim(start,final) 
    q1.grid();q2.grid();q3.grid() 
    q1.tick_params(axis='x', labelsize=16)
    q1.tick_params(axis='y', labelsize=16)
    q2.tick_params(axis='x', labelsize=16)
    q2.tick_params(axis='y', labelsize=16)
    q3.tick_params(axis='y', labelsize=16)
    q3.set_xlabel('Distance (km)')
    fig2.savefig(path+'frame='+str(i)+'_fig2.png')
cc=-1;ff1=[]
for rr,oo in enumerate(w2):
    if cc*oo<0:
        ff1.append(ox[rr])
    cc = oo
mm=-1;ff2=[]
for pp,uu in enumerate(w3):
    if mm*uu<0:
        ff2.append(ox[pp])
    mm = uu  
if len(ff2)>1 and (ff2[1]-ff2[0])>100 and ff2[0]>start:
    find_flat_dz2.append(fl.time[i])
    if len(ff1)>1 and (ff1[1]-start)>50:
        find_flat_dz1.append(i)
    
# filename2='/home/jiching/geoflac/data/'+model+'_flat_slab_time2'
# f = open(filename2 ,'w')
# for trep in range(len(find_flat_dz2)):
#     f.write('%f\n'%find_flat_dz2[trep])
# f.close()
#filename1='/home/jiching/geoflac/data/'+model+'_flat_slab_time1'
#f = open(filename1 ,'w')
#for trep in range(len(find_flat_dz1)):
#    f.write('%f\n'%find_flat_dz1[trep])
#f.close()
print(find_flat_dz2)
print(find_flat_dz1)