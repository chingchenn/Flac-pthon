#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

@author: ji-chingchen
"""


import flac
import math
import sys, os
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import cm
from heapq import nsmallest
from scipy import interpolate
import function_for_flac as fd
import function_savedata as fs
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.figsize"] = (10,12)
model = sys.argv[1]
# model = 'Ref_Nazca'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
# path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
# savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
# figpath = '/Users/chingchen/Desktop/figure/'

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
g = 10

import time
start_time = time.time()


def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt
###----------------------- Slab sinking force with time-------------------------------
#    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 
def slab_sinking_force(frame):
    # ------ read data from model -----
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    density = fl.read_density(frame)
    area = fl.read_area(frame)
    temp = fl.read_temperature(frame)
    ele_temp = temp_elements(temp)
    # ----- empty array and data -----
    rho_diff = np.zeros(nez)
    slab_area = np.zeros(nez)
    Fsb = 0 
    
    for z_ind in range(1,len(ele_z[0])):
        ind_eclogite = (phase[:,z_ind] == phase_eclogite) + (phase[:,z_ind] == phase_eclogite_1)\
            #+ (phase[:,z_ind] == phase_hydratedmantle) + (phase[:,z_ind] == phase_mantle2) # denser plate 
        man_eclogite = (phase[:,z_ind] == phase_mantle1) + (phase[:,z_ind] == phase_serpentinite)\
            + (phase[:,z_ind] == phase_hydratedmantle) # lighter plate
        ind_slab = (ele_temp[:,z_ind]<700)*((phase[:,z_ind]== phase_mantle1)+ (phase[:,z_ind]==phase_eclogite)\
                + (phase[:,z_ind] == phase_eclogite_1)+ (phase[:,z_ind] == phase_hydratedmantle))
        if True in ind_eclogite and True in man_eclogite:
            den_mantle = np.average(density[man_eclogite,z_ind])
            den_eco = np.average(density[ind_eclogite,z_ind])
            if den_eco < den_mantle:
                continue
            rho_diff[z_ind] = den_eco - den_mantle
            slab_area[z_ind] = area[ind_slab,z_ind].sum()
            #if True in ind_slab:
                #print(ele_x[ind_slab],z_ind)
        Fsb += rho_diff[z_ind] * g * slab_area[z_ind]
    return Fsb # N/m (2D)
def slab_sinking_torque(frame):
    # ------ read data from model -----
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    temp = fl.read_temperature(frame)
    temp_ele = temp_elements(temp)
    phase = fl.read_phase(frame)
    density = fl.read_density(frame)
    area = fl.read_area(frame)
    ind_trench = int(trench_index[frame])
    moment_point_x,moment_point_z  = trench_x[frame], trench_z[frame]
    ref_den = density[-4,:]
    # ----- empty array and data -----
    rho_diff = np.zeros(nex)
    Fsb = 0 
    # ----- Start Calculation ------
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
        # Choose the eclogite area
        ind_eclogite = (phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic)
        if True in ind_eclogite :
            for ele_index in range(len(ele_x[x_ind,:][ind_eclogite])):
                x1 = ele_x[x_ind,:][ind_eclogite][ele_index]
                z1 = ele_z[x_ind,:][ind_eclogite][ele_index]
                torque_length= (x1-moment_point_x) *1e3
                filter_man = (abs(ele_z[x_ind,:]-z1)<40)
                if not True in filter_man:
                    continue
                volume = np.sum(area[x_ind,:][temp_ele[x_ind,:]<800])
                rho_diff[x_ind] = np.average(density[x_ind,:]-ref_den)
                # rho_diff[x_ind] = np.average(density[x_ind,:][filter_man]-ref_den[filter_man])
                Fsb+= torque_length*rho_diff[x_ind]*g*volume
    return Fsb # N (2D)

###---------------------- Mantle flow traction force with time -------------------------------
#   Ft = sum(shear stress * dx) = \sigma_xz * dx 
def mantle_traction_force(frame):
    x,z, = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x,z)
    sxz = fl.read_sxz(frame)
    phase = fl.read_phase(frame)
    dx = np.zeros(nex)
    stressxz = np.zeros(nex)
    Ft = 0
    ## Base of the upper plate
    for qq in range(1,nex):
        upper_plate = (phase[qq,:]== phase_uppercrust) + (phase[qq,:] == phase_lowercrust)
        if True in upper_plate:
            last_deep= np.argmin(ele_z[qq,upper_plate])
            dx[qq] = (x[qq+1,last_deep]-x[qq,last_deep])*1e3 # km --> m
            stressxz[qq] = abs(sxz[qq,last_deep]*1e8) # stress kb --> N/m^2
        Ft += stressxz[qq] * dx[qq] # N/m
    return Ft # N/m (2D)

    
def suction_force(frame):    
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    pressure = fl.read_pres(frame)
    area = fl.read_area(frame)
    
    aatop = 0;aasub = 0
    Fsub = 0; Ftop = 0
    onepre=-pressure.flatten()
    a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(onepre-fit).reshape(len(phase),len(phase[0]))*1e8  # N/m^2
    final_ind = int(trench_index[frame])
    ind_trench = int(trench_index[frame])
    xoceanic = np.zeros(len(ele_z)-ind_trench)
    iiind = np.zeros(len(ele_z)-ind_trench)
    length = 0
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
        ind_oceanic = (phase[x_ind,:] == phase_oceanic) + (phase[x_ind,:] == phase_eclogite)+(phase[x_ind,:] == phase_oceanic_1) + (phase[x_ind,:] == phase_eclogite_1)
        subducted_sed = (phase[x_ind,:] == phase_sediment) + (phase[x_ind,:] ==phase_sediment_1) + (phase[x_ind,:] == phase_schist)
    
        if True in (ind_oceanic+subducted_sed):
            if (x_ind - final_ind) > 2:
                break
            final_ind = x_ind
    
            oceanic_plate_index = [ww for ww, x in enumerate(ind_oceanic+subducted_sed) if x]
            xoceanic[ii] = int(np.median(oceanic_plate_index))
            xstd = np.std(oceanic_plate_index)
            oo = oceanic_plate_index
            if xstd > 15:
                oo = np.array(oceanic_plate_index)[abs(ele_z[x_ind,oceanic_plate_index]-ele_z[x_ind,int(xoceanic[ii-1])])<30]
                if len(oo)==0:
                    continue
                xoceanic[ii] = int(np.median(oo))
            av_oc_ind = int(np.median(oo))
            iiind[ii] = av_oc_ind 
            # make sure the top element is continent
            if phase[x_ind,0]!=2 and phase[x_ind,0]!=14 and phase[x_ind,0]!=6: 
                continue
            # area of sub mantle 
            submantle = (ele_z[x_ind,av_oc_ind:]> -660)*((phase[x_ind,av_oc_ind:]==phase_mantle1) + \
                (phase[x_ind,av_oc_ind:]==phase_mantle2) + \
                (phase[x_ind,av_oc_ind:]==phase_serpentinite)+\
                (phase[x_ind,av_oc_ind:]==phase_hydratedmantle))
    
            # area of top mantle
            topmantle = (phase[x_ind,:av_oc_ind]==phase_mantle1) + \
                (phase[x_ind,:av_oc_ind]==phase_mantle2) + \
                (phase[x_ind,:av_oc_ind]==phase_serpentinite)
        
            if True in topmantle or True in submantle:
                aasub += (area[x_ind,av_oc_ind:][submantle]).sum()
                aatop += (area[x_ind,:av_oc_ind][topmantle]).sum()

                Fsub += ((dypre*area)[x_ind,av_oc_ind:][submantle]).sum()
                Ftop += ((dypre*area)[x_ind,:av_oc_ind][topmantle]).sum()                
                if ii > 1:
                    dx = ele_x[x_ind,av_oc_ind]-ele_x[x_ind-1,int(iiind[ii-1])]
                    dz = ele_z[x_ind,av_oc_ind]-ele_z[x_ind-1,int(iiind[ii-1])]
                    length += np.sqrt(dx**2 + dz**2)*1000
    if length != 0 and aasub !=0 and aatop !=0:
        Fsu = (Fsub/aasub-Ftop/aatop)*length
    else:
        Fsu = 0
    return Fsu # N/m(2D)

def find_slab_median_index2(i):    
    bet = 2
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
        
        if ox[yy-1]< 730:
            oz[yy] = np.min(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
            
        else:
            # print(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
            oz[yy] = np.max(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        px = xx

    oxx=ox[ox>start]
    oz=oz[ox>start]
    ox=oxx
    return ox,oz

def find_slab_median_index(frame):
    final_ind = int(trench_index[frame])
    ind_trench = int(trench_index[frame])
    xoceanic = np.zeros(len(ele_z)-ind_trench)
    slab_x = np.zeros(len(ele_z)-ind_trench)
    slab_z = np.zeros(len(ele_z)-ind_trench)
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
    # for ii,x_ind in enumerate(range(140,272)):
        ind_oceanic = (phase[x_ind,:] == phase_oceanic) + (phase[x_ind,:] == phase_eclogite)+(phase[x_ind,:] == phase_oceanic_1) + (phase[x_ind,:] == phase_eclogite_1)
        subducted_sed = (phase[x_ind,:] == phase_sediment) + (phase[x_ind,:] ==phase_sediment_1) + (phase[x_ind,:] == phase_schist)
        if True in (ind_oceanic+subducted_sed):
            if (x_ind - final_ind) > 2:
                break
            final_ind = x_ind
    
            oceanic_plate_index = [ww for ww, x in enumerate(ind_oceanic+subducted_sed) if x]
            xoceanic[ii] = int(np.median(oceanic_plate_index))
            xstd = np.std(oceanic_plate_index)
            oo = oceanic_plate_index
            if xstd > 15:
                print(x_ind)
                oo = np.array(oceanic_plate_index)[oceanic_plate_index>=xoceanic[ii-1]]
                if len(oo)==0:
                    continue
                xoceanic[ii] = int(np.median(oo))
            av_oc_ind = int(np.median(oo))
            slab_x[ii] = ele_x[x_ind,av_oc_ind]
            slab_z[ii] = ele_z[x_ind,av_oc_ind]
        # print(ii,av_oc_ind,oo,x_ind)
    return slab_x,slab_z

def dynamics_pressure(frame):
    pre = -fl.read_pres(frame) *1e8
    ooone = pre.flatten()
    x,z = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x, z)
    a,b=np.polyfit(pre[ele_z<-50],ele_z[ele_z<-50].flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(ooone-fit).reshape(len(pre),len(pre[0])) 
    return x,z,dypre

###---------------------- Mantle suction force with time -------------------------------
def suction_force3(frame):
    depth1 = -350
    depth2 = -10
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    x,z,dpre = dynamics_pressure(frame) # N/m^2
    slab_x,slab_z = find_slab_median_index2(frame)

    xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    
    z5 =  np.polyfit(xslab,zslab,5)
    p5 = np.poly1d(z5)
    w5 = p5(xslab) # 5st poly
    zslab = w5
    list_subslab =  [[] for i in range(len(zslab))]
    list_topslab =  [[] for i in range(len(zslab))]
    ind_trench = int(trench_index[frame])
    
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
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
        #Choose the top mantle area
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
    
    Ptotal = np.zeros(len(xslab))
    dl = np.zeros(len(xslab))
    length = np.zeros(len(xslab))
    for ss in range(len(xslab)-1):
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
            if np.linalg.norm(vec)==0:
                print(ss,x1,z1)
                continue
            cosan = np.dot(vec,ref)/np.linalg.norm(vec)
            length[ss] = dl[ss]*cosan
        if len(list_subslab[ss])!=0:
            pres = np.zeros(len(list_subslab[ss]))
            for rr in range(len(list_subslab[ss])):
                pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
                psub = np.average(pres)  
            
        if len(list_topslab[ss])!=0:
            pret = np.zeros(len(list_topslab[ss]))
            for rr in range(len(list_topslab[ss])):
                pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
                ptop = np.average(pret)
        Ptotal[ss] = psub-ptop
    Fz = (Ptotal*length).sum() # N/m
    print('Fz=',Fz/1e12)
    return Fz

def suction_force3_torque(frame):
    depth1 = -150
    depth2 = -10
    dis_range = 20
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    x,z,dpre = dynamics_pressure(frame) # N/m^2
    slab_x,slab_z = find_slab_median_index2(frame)
    
    xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    
    z5 =  np.polyfit(xslab,zslab,5)
    p5 = np.poly1d(z5)
    w5 = p5(xslab) # 5st poly
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

            for ele_index in range(len(ele_x[x_ind,:][submantle])):
                maxmax = 9999
                x1 = ele_x[x_ind,:][submantle][ele_index]
                z1 = ele_z[x_ind,:][submantle][ele_index]              
                if not True in (np.sqrt((xslab-x1)**2+(zslab-z1)**2)<dis_range+1):
                    continue
                for slab_ind in range(len(xslab)):
                    dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                    if dis > dis_range:
                        continue
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
                    if dis > dis_range:
                        continue
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
    
    Ptotal = np.zeros(len(xslab))
    dl = np.zeros(len(xslab))
    torque_length = np.zeros(len(xslab))
    P0 = np.array((trench_x[frame], trench_z[frame]))
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
                pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
                psub = np.average(pres)
        if len(list_topslab[ss])!=0:
            pret = np.zeros(len(list_topslab[ss]))
            for rr in range(len(list_topslab[ss])):
                pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
                ptop = np.average(pret)
        Ptotal[ss] = psub-ptop
        Fsu += Ptotal[ss]*dl[ss]*torque_length[ss]*costheta
        # Fz = (Ptotal*length).sum() # N/m
        # print('Fz=',Fz/1e12)
        # print('Fsu=',Fsu/1e12)
    return Fsu

###---------------------- couple zone -------------------------------
def shearstress_indistance(frame):
    phase_uppercrust = 2
    phase_lowercrust = 14
    x,z, = fl.read_mesh(frame)
    ele_x,ele_z = flac.elem_coord(x,z)
    phase = fl.read_phase(frame)
    sxz = fl.read_sxz(frame)
    stressxz = np.zeros(nex)
    for qq in range(1,nex):
        upper_plate = (phase[qq,:]== phase_uppercrust) + (phase[qq,:] == phase_lowercrust)
        if True in upper_plate:
            last_deep= np.argmin(ele_z[qq,upper_plate])
            stressxz[qq] = sxz[qq,last_deep]*1e2
    ssxz = fd.moving_window_smooth(stressxz,8)
    return ele_x[:,0], ssxz

# fig, (ax2) = plt.subplots(1,1,figsize=(15,9))
# rainbow = cm.get_cmap('gray_r',end)
# newcolors = rainbow(np.linspace(0, 1, end))
# for i in range(40,end,20):
#     dis, ssxz = shearstress_indistance(i)
#     ax2.plot(dis-trench_x[i],ssxz,c = newcolors[i],lw=4,label=str(round(fl.time[i],1))+' Myr')
# ax2.set_xlim(0,700)
# #ax2.set_xlim(np.average(trench_x),1200)
# ax2.set_title(model+' $\sigma_{xz}$ and distance',fontsize = 24)
# ax2.set_xlabel('Distance from trench (km)',fontsize=16)
# ax2.set_ylabel('$\sigma_{xz}$ (MPa)',fontsize=16)
# ax2.tick_params(axis='x', labelsize=16)
# ax2.tick_params(axis='y', labelsize=16)
# ax2.grid()
# ax2.legend(fontsize=20)
# 
# ax2.spines['bottom'].set_linewidth(bwith)
# ax2.spines['top'].set_linewidth(bwith)
# ax2.spines['right'].set_linewidth(bwith)
# ax2.spines['left'].set_linewidth(bwith)

# #fig.savefig('/home/jiching/geoflac/figure/'+model+'_sxz_dis.png')

#------------------------------------------------------------------------------
if __name__ == '__main__':
# if __name__ == '__kk__':
    fsb=np.zeros(end)
    ft=np.zeros(end)
    fsu=np.zeros(end)
    ratio=np.zeros(end)
    
    for i in range(5,end):
        loop_time = time.time()
        # fsb[i] = slab_sinking_force(i)
        fsb[i] = slab_sinking_torque(i)
        fsu[i] = suction_force3_torque(i)

    fs.save_3txt(model+'_forces',savepath,fl.time,fsb,fsu)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    fl.time,fbb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    
    sb = fd.moving_window_smooth(fsb,8)
    tt = fd.moving_window_smooth(fsu,8)
    ax.plot(fl.time,sb,c='#c06c84',label='slab pull (N)',lw=4)
    ax.plot(fl.time,tt,c="#355c7d",label='suction force (N)',lw=4)
    ax.legend(fontsize=16,loc='upper left')
    #================================figure setting================================
    ax.set_xlabel('Time (Myr)',fontsize=16)
    ax.set_ylabel('Force (N/m)',fontsize=16)
    ax.set_xlim(0, fl.time[-1])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid()
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    #ax.set_yscale('log')
    ax.set_title('Forces of '+model,fontsize=20)
    fig.savefig(figpath+model+'_slab_torque.png')
    # fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
    # ratio_f = fd.moving_window_smooth(ratio[ratio>0],5)
    # ax2.plot(fl.time[ratio>0],ratio_f,c="#355c7d",label='ratio of these forces)',lw=4)
    # ax2.set_xlabel('Time (Myr)',fontsize=16)
    # ax2.set_xlim(0, fl.time[-1])
    # ax2.set_ylim(0, 30)
    # ax2.tick_params(axis='x', labelsize=16)
    # ax2.tick_params(axis='y', labelsize=16)
    # ax2.grid()
    # ax2.spines['bottom'].set_linewidth(bwith)
    # ax2.spines['top'].set_linewidth(bwith)
    # ax2.spines['right'].set_linewidth(bwith)
    # ax2.spines['left'].set_linewidth(bwith)
    # fig2.savefig(figpath+model+'_slab_force_ratio.png')
print("--- %s seconds ---" % (time.time() - start_time)/60)
