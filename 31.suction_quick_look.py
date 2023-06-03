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
model = 'Nazca_a0702'
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
x_limit = 1000
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

fig2 = 0
fig3 = 1


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


def find_mid_point(x,z):
    midx=np.zeros(len(x)-1)
    midz=np.zeros(len(z)-1)
    for qq in range(len(x)-1):
        midx[qq] = 0.5*(x[qq]+x[qq+1])
        midz[qq] = 0.5*(z[qq]+z[qq+1])
    return midx,midz


def find_slab_median_index2(i):    
    bet = 0.5
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
frame = 180
dis_range = 20
depth1 = -150
depth2 = -10

bwith = 3
for frame in [30,90,150]:
# for frame in [180]:

    ###---------------------------------------------------------------------------------------------------------
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    _,_,dpre = dynamics_pressure(frame) # N/m^2
    slab_x,slab_z = find_slab_median_index2(frame)
    density = fl.read_density(frame)
    
    
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
    # fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    # # ax.pcolormesh(ele_x,-ele_z,dpre/1e6,cmap=cm,vmin=-200, vmax=200)
    # # ax.pcolormesh(x,-z,density,cmap=cm,vmin=2900,vmax=3500)
    # ax.pcolormesh(ele_x,-ele_z,phase,cmap=phase15,vmin=1, vmax=20)
    # ax.set_ylim(300,-0)
    # ax.set_xlim(200,x_limit)
    # # ax.set_xlim(trench_x[frame]-200,min(trench_x[frame]+800,1200))
    # ax.set_aspect('equal')
    # # xt,zt = fl.read_mesh(frame)
    # temp = fl.read_temperature(frame)
    # ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)

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
    ind_trench = int(trench_index[frame])
    moment_point_x,moment_point_z  = trench_x[frame], trench_z[frame]
    # ax.scatter(moment_point_x,-moment_point_z,c = 'green')
    ref_den = density[-4,:]
    # ----- empty array and data -----
    Fsb = 0 
    Fsbx = np.zeros(len(ele_x))
    force_sb = np.zeros(len(ele_x))
    
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
    # for ii,x_ind in enumerate(range(166,182)):
        # Choose the eclogite area
        ind_eclogite = (ele_z[x_ind,:]<-20)*((phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic))
        if not True in ind_eclogite:
            # print(frame,x_ind)
            continue
        top_slab_index = np.where(ind_eclogite)[0][0]
        litho800 = (temp_ele[x_ind,:]<800)*(ele_z[x_ind,:]<ele_z[x_ind,top_slab_index])
        if True in litho800:
            # ax.scatter(ele_x[x_ind,:][litho800],-ele_z[x_ind,:][litho800],c = 'yellow',s=4)
            for ele_index in range(len(ele_x[x_ind,:][litho800])):
                x1 = ele_x[x_ind,:][litho800][ele_index]
                # z1 = ele_z[x_ind,:][litho800][ele_index]
                rho_diff = density[x_ind,:][litho800][ele_index]-ref_den[ele_index]
                torque_length= (x1-moment_point_x) *1e3        
                volume = area[x_ind,:][litho800][ele_index]
                Fsb+= torque_length*rho_diff*g*volume
                Fsbx[x_ind] = Fsb
                force_sb[x_ind] = rho_diff*g*volume
        # return Fsb # N (2D)
    #---------------------------------------------------------------------------
    trx = trench_x[frame]
    zslab = w5[xslab>trx]
    xslab = xslab[xslab>trx]
    # ax.scatter(xslab,-zslab,c = 'k',s = 30)
    midslabx,midslabz = find_mid_point(xslab,zslab)
    list_subslab =  [[] for i in range(len(midslabz))]
    list_topslab =  [[] for i in range(len(midslabz))]
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
                if not True in (np.sqrt((midslabx-x1)**2+(midslabz-z1)**2)<dis_range+1):
                    continue
                for slab_ind in range(len(midslabx)):
                    dis = np.sqrt((midslabx[slab_ind]-x1)**2+(midslabz[slab_ind]-z1)**2)
                    # print(slab_ind, dis,ele_index, np.sqrt((xslab-x1)**2+(zslab-z1)**2))
                    # qq = (np.sqrt((xslab-x1)**2+(zslab-z1)**2)<11)
                    if dis > dis_range:
                        continue
                    if dis<maxmax:
                        maxmax = dis
                        choosex = midslabx[slab_ind]
                        choosez = midslabz[slab_ind]
                        chooseind = slab_ind
                if maxmax < dis_range:
                    if chooseind  == len(midslabx)-1:
                        kk = len(midslabx)-1
                    else:
                        kk = chooseind+1
                    pointp = np.array([x1,z1])
                    pointa = np.array([choosex,choosez])
                    pointb = np.array([midslabx[kk],midslabz[kk]])
                    if isabove(pointp,pointa,pointb)==False:
                        z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                        list_subslab[chooseind].append([x_ind,z_ind])
                        # ax.scatter(ele_x[x_ind,z_ind],-ele_z[x_ind,z_ind],c='darkred',s = 5)
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
                if not True in (np.sqrt((midslabx-x1)**2+(midslabz-z1)**2)<dis_range+1):
                    continue
                for slab_ind in range(len(midslabx)):
                    dis = np.sqrt((midslabx[slab_ind]-x1)**2+(midslabz[slab_ind]-z1)**2)
                    if dis > dis_range:
                        continue
                    if dis<maxmax:
                        maxmax = dis
                        choosex = midslabx[slab_ind]
                        choosez = midslabz[slab_ind]
                        chooseind = slab_ind
                if maxmax < dis_range:
                    if chooseind  == len(midslabx)-1:
                        kk = len(midslabx)-1
                    else:
                        kk = chooseind+1
                    pointp = np.array([x1,z1])
                    pointa = np.array([choosex,choosez])
                    pointb = np.array([midslabx[kk],midslabz[kk]])
                    if isabove(pointp,pointa,pointb)==True:
                        z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                        list_topslab[chooseind].append([x_ind,z_ind])
                        # ax.scatter(ele_x[x_ind,z_ind],-ele_z[x_ind,z_ind],c='#00FF00',s = 5) 
    
    Ptotal = np.zeros(len(midslabx))
    dl = np.zeros(len(midslabx))
    torque_length = np.zeros(len(midslabx))
    P0 = np.array((trench_x[frame], trench_z[frame]))*1000
    Fsux2 = np.zeros(len(xslab))
    beta = np.zeros(len(xslab))
    force_su = np.zeros(len(xslab))
    for ss in range(len(midslabx)):
    # for ss in range(263,264):    
        psub = 0
        ptop = 0
        costheta= 0
        if ss >= 1:
            P3 = np.array((midslabx[ss-1],midslabz[ss-1]))*1e3 # midpoint
            P2 = np.array((xslab[ss],zslab[ss]))*1e3 # 
            P1 = np.array((xslab[ss-1],zslab[ss-1]))*1e3
            dl[ss] = np.linalg.norm(P2-P1)
            torque_length[ss]= np.linalg.norm(P3-P0)
            PPP = np.array([P3-P0,P3-P1])
            costheta = np.dot(PPP[0],PPP[1])/(dl[ss]/2)/torque_length[ss]
            beta[ss] = np.arccos(costheta)*180/np.pi
        if len(list_subslab[ss])!=0:
            pres = np.zeros(len(list_subslab[ss]))
            for rr in range(len(list_subslab[ss])):
                indx = list_subslab[ss][rr][0]
                indz = list_subslab[ss][rr][1]
                xm = ele_x[indx,indz]
                zm = ele_z[indx,indz]
                # ax.scatter(xm,-zm,c='b',s = 10)
                pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
                psub = np.average(pres)
        if len(list_topslab[ss])!=0:
            pret = np.zeros(len(list_topslab[ss]))
            for rr in range(len(list_topslab[ss])):
                indx = list_topslab[ss][rr][0]
                indz = list_topslab[ss][rr][1]
                xm = ele_x[indx,indz]
                zm = ele_z[indx,indz]
                # ax.scatter(xm,-zm,c='r',s = 10)
                pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
                ptop = np.average(pret)
        Ptotal[ss] = psub-ptop
        Fsu += Ptotal[ss]*dl[ss]*torque_length[ss]*costheta
        Fsux2[ss] = Fsu
        force_su[ss] = Ptotal[ss]*dl[ss]
    if fig2:
        fig2, (ax2,ax3,ax4)= plt.subplots(3,1,figsize=(15,10))
        ax2.plot(ele_x[:,0][Fsbx!=0],Fsbx[Fsbx!=0],c ='#524B52',label='slab pull (N)',lw=4) 
        ax2.plot(xslab[Fsux2!=0],Fsux2[Fsux2!=0],c="#DC143C",label='suction torque (N)',lw=4)
        ax3.scatter(xslab,beta,c='darkgreen')
        ax4.scatter(ele_x[:,0][force_sb!=0],force_sb[force_sb!=0],c ='#524B52',label='slab pull force')
        ax4.scatter(xslab[force_su!=0],force_su[force_su!=0],c="#DC143C",label = 'suction force')
        for axx in [ax2,ax3,ax4]:
            axx.grid()
            axx.set_xlim(200,x_limit)
            #axx.set_title(model+' '+str(round(frame*0.2))+' Myr',fontsize=24)
            axx.spines['bottom'].set_linewidth(bwith)
            axx.spines['top'].set_linewidth(bwith)
            axx.spines['right'].set_linewidth(bwith)
            axx.spines['left'].set_linewidth(bwith)
            axx.tick_params(axis='x', labelsize=12)
            axx.tick_params(axis='y', labelsize=12)
            
        ax2.legend(fontsize=16,loc='upper left')
        ax4.legend(fontsize=16,loc='upper left')
        ax2.set_ylabel('Torque (N)',fontsize=16)
        ax3.set_ylabel('Angle (degree)',fontsize=16)
        ax4.set_ylabel('force (N/m)',fontsize=16)
        ax4.set_xlabel('X Distance',fontsize=16)


    if fig3:
        fig3, (ax,aa,ax2,ax4)= plt.subplots(4,1,figsize=(8,15))
        ax.pcolormesh(ele_x,-ele_z,phase,cmap=phase15,vmin=1, vmax=20)
        ax.set_ylim(300,-0)
        ax.set_aspect('equal')
        temp = fl.read_temperature(frame)
        ax.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)
        aa.pcolormesh(ele_x,-ele_z,dpre/1e6,cmap=cm,vmin=-200, vmax=200)
        aa.set_ylim(300,-0)
        aa.set_aspect('equal')
        temp = fl.read_temperature(frame)
        aa.contour(x,-z,temp,colors='0.5',levels =[200,400,600,800,1000,1200],linewidths=3)

        ax2.plot(ele_x[:,0][Fsbx!=0],Fsbx[Fsbx!=0],c ='#524B52',label='slab pull (N)',lw=4) 
        ax2.plot(xslab[Fsux2!=0],Fsux2[Fsux2!=0],c="#DC143C",label='suction torque (N)',lw=4)
        ax4.scatter(ele_x[:,0][force_sb!=0],force_sb[force_sb!=0],c ='#524B52',label='slab pull force')
        ax4.scatter(xslab[force_su!=0],force_su[force_su!=0],c="#DC143C",label = 'suction force')
        for axx in [ax,aa,ax2,ax4]:
            axx.set_xlim(200,x_limit)
            #axx.set_title(model+' '+str(round(frame*0.2))+' Myr',fontsize=24)
            axx.spines['bottom'].set_linewidth(bwith)
            axx.spines['top'].set_linewidth(bwith)
            axx.spines['right'].set_linewidth(bwith)
            axx.spines['left'].set_linewidth(bwith)
            axx.tick_params(axis='x', labelsize=12)
            axx.tick_params(axis='y', labelsize=12)
        for aaa in [ax2,ax4]:
            aaa.grid()
            
        ax2.legend(fontsize=16,loc='upper left')
        ax4.legend(fontsize=16,loc='lower left')
        ax2.set_ylabel('Torque (N)',fontsize=16)
        ax4.set_ylabel('force (N/m)',fontsize=16)
        ax4.set_xlabel('X Distance',fontsize=16)