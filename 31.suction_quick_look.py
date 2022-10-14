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
from heapq import nlargest,nsmallest
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.interpolate import  UnivariateSpline,Akima1DInterpolator, PchipInterpolator
#------------------------------------------------------------------------------
# plt.rcParams["font.family"] = "Times New Roman"
# plt.rcParams["figure.figsize"] = (10,12)
#model = sys.argv[1]
model = 'ch1520'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
path = 'D:/model/'
path = '/Users/chingchen/Desktop/model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
figpath='/scratch2/jiching/figure/'
figpath = '/Users/chingchen/Desktop/figure/'
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


def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt

def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre

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

g = 10
frame = 140
dis_range = 25
depth1 = -150
depth2 = -10

bwith = 3
for frame in [120]:



    ###---------------------------------------------------------------------------------------------------------
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    _,_,dpre = dynamics_pressure(frame) # N/m^2
    slab_x,slab_z = find_slab_median_index2(frame)
    
    xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    
    z5 =  np.polyfit(xslab,zslab,5)
    p5 = np.poly1d(z5)
    w5 = p5(xslab) # 5st poly
    ###----------------------------------------------------------------------------
    colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
              "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
              "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
              "#F67280","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    cm = plt.cm.get_cmap('RdYlBu_r')
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    ax.scatter(ele_x,ele_z,c=-dpre/1e6,cmap=cm,vmin=-200, vmax=200,s=40)
    # ax.scatter(ele_x,ele_z,c=phase,cmap=phase15,vmin=1, vmax=20,s=100)
    ax.set_ylim(-300,-0)
    ax.set_xlim(trench_x[frame]-200,min(trench_x[frame]+800,1200))
    ax.set_aspect('equal')
    #plt.scatter(xslab,zslab,c = 'k',s = 50)
    # ax.scatter(xslab,w5,c = 'k',s = 30)
    
    
    
    ###----------------------- Slab sinking force with time-------------------------------
    #    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 
    # def slab_sinking_torque(frame):
        # ------ read data from model -----
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    density = fl.read_density(frame)
    area = fl.read_area(frame)
    temp = fl.read_temperature(frame)
    temp_ele = temp_elements(temp)
    # ----- empty array and data -----
    rho_diff = np.zeros(nex)
    rho_diff2 = np.zeros(nex)
    Fsb = 0 
    Fsbx = np.zeros(len(ele_x))
    ind_trench = int(trench_index[frame])
    moment_point_x,moment_point_z  = trench_x[frame], trench_z[frame]
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
    # for ii,x_ind in enumerate(range(ind_trench,ind_trench+50)):
        # Choose the eclogite area
        ind_eclogite = (phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic)
        ref_den = density[-5,:]
        if True in ind_eclogite:# and True in man_eclogite:
            # ax.scatter(ele_x[x_ind,:][man_eclogite],ele_z[x_ind,:][man_eclogite],c = 'b')
            for ele_index in range(len(ele_x[x_ind,:][ind_eclogite])):
                x1 = ele_x[x_ind,:][ind_eclogite][ele_index]
                z1 = ele_z[x_ind,:][ind_eclogite][ele_index]
                torque_length= (x1-moment_point_x) *1e3        
                filter_man = (abs(ele_z[x_ind,:]-z1)<40)
                if not True in filter_man:
                    continue
                volume = np.sum(area[x_ind,:][temp_ele[x_ind,:]<800])
                # rho_diff[x_ind] = np.average(density[x_ind,:][filter_man]-ref_den[filter_man])
                rho_diff2[x_ind] = np.average(density[x_ind,:]-ref_den)
                Fsb+= torque_length*rho_diff2[x_ind]*g*volume
                Fsbx[x_ind] = Fsb
        # return Fsb # N (2D)
    
    # ###------------------------------------------------------------------------------------------------
    
    zslab = w5
    list_subslab =  [[] for i in range(len(zslab))]
    list_topslab =  [[] for i in range(len(zslab))]
    ind_trench = int(trench_index[frame])
    Fsu = 0
    
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
        # Choose the submantle area
        isabove = lambda p,a,b : np.cross(p-a,b-a)<0
        submantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*((phase[x_ind,:]==phase_mantle1) + \
                (phase[x_ind,:]==phase_mantle2) +(phase[x_ind,:]==phase_hydratedmantle))
        if True in submantle:
            # ax.scatter(ele_x[x_ind,:][submantle],ele_z[x_ind,:][submantle],c = 'b')
            # print(len(ele_x[x_ind,:][submantle]))
            for ele_index in range(len(ele_x[x_ind,:][submantle])):
                maxmax = 9999
                x1 = ele_x[x_ind,:][submantle][ele_index]
                z1 = ele_z[x_ind,:][submantle][ele_index]
                
                # print(np.sqrt((xslab-x1)**2+(zslab-z1)**2))
                # print(<11)
                if not True in (np.sqrt((xslab-x1)**2+(zslab-z1)**2)<dis_range+1):
                    continue
                
                for slab_ind in range(len(xslab)):
                    dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                    # print(slab_ind, dis,ele_index, np.sqrt((xslab-x1)**2+(zslab-z1)**2))
                    # qq = (np.sqrt((xslab-x1)**2+(zslab-z1)**2)<11)
                    if dis<maxmax:
                        maxmax = dis
                        choosex = xslab[slab_ind]
                        choosez = zslab[slab_ind]
                        chooseind = slab_ind
                if maxmax < dis_range:
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
                        # plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='b',s = 5)
                    #elif isabove(pointp,pointa,pointb)==True:
                        #z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                        #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 
        #Choose the top mantle area
        topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*\
        ((phase[x_ind,:]==phase_mantle1) + (phase[x_ind,:]==phase_mantle2) + \
                (phase[x_ind,:]==phase_serpentinite))
        if True in topmantle:
            for ele_index in range(len(ele_x[x_ind,:][topmantle])):
                maxmax = 9999
                x1 = ele_x[x_ind,:][topmantle][ele_index]
                z1 = ele_z[x_ind,:][topmantle][ele_index]
                if not True in (np.sqrt((xslab-x1)**2+(zslab-z1)**2)<dis_range+1):
                    continue
                for slab_ind in range(len(xslab)):
    
                    dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                    if dis<maxmax:
                        maxmax = dis
                        choosex = xslab[slab_ind]
                        choosez = zslab[slab_ind]
                        chooseind = slab_ind
                if maxmax < dis_range:
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
                        # ax.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 
    Ptotal = np.zeros(len(xslab))
    dl = np.zeros(len(xslab))
    torque_length = np.zeros(len(xslab))
    P0 = np.array((trench_x[frame], trench_z[frame]))
    Fsux = np.zeros(len(xslab))
    for ss in range(len(xslab)-1):
        psub = 0
        ptop = 0
        costheta= 0
        if ss >= 1:
            P1 = np.array((xslab[ss],zslab[ss]))*1e3
            P2 = np.array((xslab[ss-1],zslab[ss-1]))*1e3
            dl[ss] = np.linalg.norm(P1-P2)
            torque_length[ss]= np.linalg.norm(P1-P0)
            PPP = np.array([P1-P0,P1-P2])
            costheta = np.dot(PPP[0],PPP[1])/dl[ss]/torque_length[ss]
        if len(list_subslab[ss])!=0:
            pres = np.zeros(len(list_subslab[ss]))
            for rr in range(len(list_subslab[ss])):
                indx = list_subslab[ss][rr][0]
                indz = list_subslab[ss][rr][1]
                xm = ele_x[indx,indz]
                zm = ele_z[indx,indz]
                # plt.scatter(xm,zm,c='b',s = 5)
                pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
                psub = np.average(pres)
            
        if len(list_topslab[ss])!=0:
            pret = np.zeros(len(list_topslab[ss]))
            for rr in range(len(list_topslab[ss])):
                indx = list_topslab[ss][rr][0]
                indz = list_topslab[ss][rr][1]
                xm = ele_x[indx,indz]
                zm = ele_z[indx,indz]
                # plt.scatter(xm,zm,c='r',s = 5)
                pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
                ptop = np.average(pret)
        Ptotal[ss] = psub-ptop
        Fsu += Ptotal[ss]*dl[ss]*torque_length[ss]*costheta
        Fsux[ss] = Fsu
        # Fz = (Ptotal*length).sum() # N/m
        # print('Fz=',Fz/1e12)
        # print('Fsu=',Fsu/1e12)
    fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
    ax2.plot(ele_x[:,0][Fsbx!=0],Fsbx[Fsbx!=0],c ='#c06c84',label='slab pull (N)',lw=4) 
    ax2.plot(xslab[Fsux!=0],Fsux[Fsux!=0],c="#355c7d",label='suction force (N)',lw=4)
    ax2.legend(fontsize=16,loc='upper left')
    ax2.grid()
    ax2.set_title(model+' '+str(round(frame*0.2))+' Myr',fontsize=24)
    