#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:29:55 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import math
import matplotlib
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.interpolate import  UnivariateSpline,Akima1DInterpolator, PchipInterpolator
#------------------------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.figsize"] = (10,12)
#model = sys.argv[1]
model = 'b0801c'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
#path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/Users/ji-chingchen/Desktop/figure/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
nex = fl.nx-1; nez=fl.nz-1
time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T

phase_uppercrust = 2
phase_oceanic = 3
phase_mantle1 = 4
phase_schist = 5
phase_mantle2 = 8
phase_serpentinite = 9
phase_sediment = 10
phase_sediment_1 = 11
phase_eclogite = 13
phase_lowercrust = 14
phase_hydratedmantle = 16
phase_oceanic_1 = 17
phase_eclogite_1 = 18

bwith = 3
def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt
def dynamic_pressure(frame):
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    pressure = -fl.read_pres(frame)
    onepre=pressure.flatten()
    a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(onepre-fit).reshape(len(ele_x),len(ele_x[0]))*1e8  # N/m^2
    return dypre # N/m^2

def find_slab_median_index2(i):    
    bet = 1
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
    ## In this code, we considered the marker phase, not the element phase
    x_trench = trench_x[i]
    x_ocean = mx[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]
    z_ocean = mz[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]
    # if z_trench> -2 or min(z_ocean)>-200:
        # continue
    start = math.floor(x_trench-50)
    final = math.floor(np.max(x_ocean))
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    #find initial basalt depth to remove the weage basalt
    if len(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])==0:
            print('no')
    kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
    x_ocean = x_ocean[z_ocean<kk]
    z_ocean = z_ocean[z_ocean<kk]
    # interplate to the grid length "bet"
    for yy,xx in enumerate(x_grid):
        if len(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])==0:
            continue

        oz[yy] = np.min(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        px = xx
    oxx=ox[ox>start]
    oz=oz[ox>start]
    ox=oxx
    return ox,oz

depth1 = -350
depth2 = -10
frame = 140
x, z = fl.read_mesh(frame)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(frame)
dpre = dynamic_pressure(frame) # N/m^2
slab_x,slab_z = find_slab_median_index2(frame)

xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]


z1=np.polyfit(xslab,zslab,4)
p4=np.poly1d(z1)
w1=p4(xslab) # 4st poly 
p3=np.polyder(p4,1)
p2=np.polyder(p4,2)
w2=p3(xslab) # 3st poly 
w3=p2(xslab) # 2st poly

z5 =  np.polyfit(xslab,zslab,5)
p5 = np.poly1d(z5)
w5 = p5(xslab) # 5st poly

list_subslab =  [[] for i in range(len(zslab))]
list_topslab =  [[] for i in range(len(zslab))]
ind_trench = int(trench_index[frame])

colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase15= matplotlib.colors.ListedColormap(colors)
cm = plt.cm.get_cmap('RdYlBu_r')
#plt.scatter(ele_x,ele_z,c=-dpre/1e6,cmap=cm,vmin=-200, vmax=200,s=40)
plt.scatter(ele_x,ele_z,c=phase,cmap=phase15,vmin=1, vmax=20,s=100)
plt.ylim(-300,-0)
plt.xlim(300,1000)
plt.axes().set_aspect('equal')
#plt.scatter(xslab,zslab,c = 'k',s = 50)
plt.scatter(xslab,w5,c = 'k',s = 30)
zslab = w5

#p = 0
for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
#for ii,x_ind in enumerate(range(ind_trench,170)):
    zind = np.where(ele_z[x_ind,:]==slab_z[ii])[0]
    if len(zind)==0:
        zind = 0
    else:
        zind = int(zind)
    # Choose the submantle area
    isabove = lambda p,a,b : np.cross(p-a,b-a)<0
    submantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*((phase[x_ind,:]==phase_mantle1) + \
            (phase[x_ind,:]==phase_mantle2) +(phase[x_ind,:]==phase_hydratedmantle))
    if True in submantle:
        for ele_index in range(len(ele_x[x_ind,:][submantle])):
            maxmax = 9999
            x1 = ele_x[x_ind,:][submantle][ele_index]
            z1 = ele_z[x_ind,:][submantle][ele_index]
            for slab_ind in range(len(xslab)):
                dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                if dis<maxmax:
                    maxmax = dis
                    choosex = xslab[slab_ind]
                    choosez = zslab[slab_ind]
                    chooseind = slab_ind
            if maxmax < 50:
                if chooseind  == len(xslab)-1:
                    kk = len(xslab)-1
                else:
                    kk = chooseind+1
                pointp = np.array([x1,z1])
                pointa = np.array([choosex,choosez])
                pointb = np.array([xslab[kk],zslab[kk]])
                if isabove(pointp,pointa,pointb)==False:
                    z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                    list_subslab[chooseind].append([x_ind,z_ind])
                    #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='b',s = 5)
                #elif isabove(pointp,pointa,pointb)==True:
                    #z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                    #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 
    #Choose the top mantle area
    if slab_z[ii]==0:
        slab_z[ii]=-410
    #topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*\
    #((phase[x_ind,:]==phase_mantle1) + (phase[x_ind,:]==phase_mantle2) + \
      #       (phase[x_ind,:]==phase_serpentinite)+(ele_x[x_ind,:]>max(slab_x)))*(ele_z[x_ind,:]>slab_z[ii])
    topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*\
    ((phase[x_ind,:]==phase_mantle1) + (phase[x_ind,:]==phase_mantle2) + \
            (phase[x_ind,:]==phase_serpentinite))
    if True in topmantle:
        for ele_index in range(len(ele_x[x_ind,:][topmantle])):
            maxmax = 9999
            x1 = ele_x[x_ind,:][topmantle][ele_index]
            z1 = ele_z[x_ind,:][topmantle][ele_index]
            for slab_ind in range(len(xslab)):
                dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                if dis<maxmax:
                    maxmax = dis
                    choosex = xslab[slab_ind]
                    choosez = zslab[slab_ind]
                    chooseind = slab_ind
            if maxmax < 50:
                if chooseind  == len(xslab)-1:
                    kk = len(xslab)-1
                else:
                    kk = chooseind+1
                pointp = np.array([x1,z1])
                pointa = np.array([choosex,choosez])
                pointb = np.array([xslab[kk],zslab[kk]])
                if isabove(pointp,pointa,pointb)==True:
                    z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                    list_topslab[chooseind].append([x_ind,z_ind])
                    #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 

Ptotal = np.zeros(len(xslab))
dl = np.zeros(len(xslab))
length = np.zeros(len(xslab))
for ss in range(len(xslab)-1):
#for ss in range(230,233):
    psub = 0
    ptop = 0
    if ss >= 1:
        x0 = xslab[ss]
        z0 = zslab[ss]
        x1 = xslab[ss-1]
        z1 = zslab[ss-1]
        dl[ss] = np.sqrt((-(zslab[ss]-zslab[ss-1]))**2+(-(xslab[ss]-xslab[ss-1]))**2)*1e3  # m
        vec = [x0-x1,z0-z1]
        ref = [1,0]
        #print(vec)
        if np.linalg.norm(vec)==0:
            print(ss,x1,z1)
            continue
        #sinan = np.cross(vec,ref)/np.linalg.norm(vec)
        cosan = np.dot(vec,ref)/np.linalg.norm(vec)
        length[ss] = dl[ss]*cosan
        
    pres = np.zeros(len(list_subslab[ss]))
    pret = np.zeros(len(list_topslab[ss]))
    for rr in range(len(list_subslab[ss])):
        indx = list_subslab[ss][rr][0]
        indz = list_subslab[ss][rr][1]
        xm = ele_x[indx,indz]
        zm = ele_z[indx,indz]
        plt.scatter(xm,zm,c='b',s = 5)
        pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
    if len(list_subslab[ss]) !=0:
        psub = np.average(pres)
    for rr in range(len(list_topslab[ss])):
        indx = list_topslab[ss][rr][0]
        indz = list_topslab[ss][rr][1]
        xm = ele_x[indx,indz]
        zm = ele_z[indx,indz]
        plt.scatter(xm,zm,c='r',s = 5)
        pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
    if len(list_topslab[ss]) !=0:
        ptop = np.average(pret)
    Ptotal[ss] = psub-ptop
Fz = (Ptotal*length).sum()
print('Fz=',Fz/1e13)