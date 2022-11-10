#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

@author: ji-chingchen
"""


import time
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
#model = 'Ref_Cocos'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
#path = '/scratch2/jiching/22summer/'
#path = '/scratch2/jiching/03model/'
#path = '/scratch2/jiching/04model/'
path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
savepath = '/Users/chingchen/Desktop/data/'
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

bwith = 3
g = 10
import time
start_time = time.time()

def find_mid_point(x,z):
    midx=np.zeros(len(x)-1)
    midz=np.zeros(len(z)-1)
    for qq in range(len(x)-1):
        midx[qq] = 0.5*(x[qq]+x[qq+1])
        midz[qq] = 0.5*(z[qq]+z[qq+1])
    return midx,midz

def temp_elements(temp):
    ttt = (temp[:fl.nx-1,:fl.nz-1] + temp[1:,:fl.nz-1] + temp[1:,1:] + temp[:fl.nx-1,1:]) / 4.
    return ttt
###----------------------- Slab sinking force with time-------------------------------
#    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 

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
    Fsbx = np.zeros(len(ele_x))
    Fsbxs = np.zeros(len(ele_x))
    # ----- Start Calculation ------
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
        # Choose the eclogite area
        ind_eclogite = (phase[x_ind,:] == phase_eclogite) + (phase[x_ind,:] == phase_eclogite_1) + (phase[x_ind,:] == phase_oceanic)
        if not True in ind_eclogite:
            continue
        top_slab_index = np.where(ind_eclogite)[0][0]
        litho800 = (temp_ele[x_ind,:]<700)*(ele_z[x_ind,:]<ele_z[x_ind,top_slab_index])
        if True in litho800 :
            for ele_index in range(len(ele_x[x_ind,:][litho800])):
                x1 = ele_x[x_ind,:][litho800][ele_index]
                rho_diff = density[x_ind,:][litho800][ele_index]-ref_den[ele_index]
                torque_length= (x1-moment_point_x) *1e3        
                volume = area[x_ind,:][litho800][ele_index]
                Fsb+= torque_length*rho_diff*g*volume
                Fsbx[x_ind] = Fsb
                Fsbxs[x_ind] = torque_length*rho_diff*g*volume
    if frame < 10:
        qq = '00'+str(frame)
    elif frame < 100 and frame >=10:
        qq = '0'+str(frame)
    else:
        qq=str(frame)
    fs.save_3txt(model+'_gravityx_'+str(qq),savepath,ele_x[:,0][Fsbx!=0],Fsbx[Fsbx!=0],Fsbxs[Fsbx!=0])
    
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

    
def find_slab_median_index2(i):    
    bet = 2
    x, z = fl.read_mesh(i)
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
    ## In this code, we considered the marker phase, not the element phase
    x_trench = trench_x[i]
    x_ocean = mx[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]
    z_ocean = mz[((phase==phase_eclogite)+(phase==phase_oceanic))*(mz>-300)]

    start = math.floor(x_trench-50)
    final = math.floor(np.max(x_ocean))
    x_grid = np.arange(start,final,bet)
    ox = np.zeros(len(x_grid))
    oz = np.zeros(len(x_grid))
    px = start-bet
    #find initial basalt depth to remove the weage basalt
    # print(i,z_ocean[(x_ocean>=start)*(x_ocean<=start+bet)])
    # print(len(z_ocean[(x_ocean>=start)*(x_ocean<=start+bet)]))
    # if len(z_ocean[(x_ocean>=start)*(x_ocean<=start+bet)])==0:
        # print('no')
    # kk=np.max(z_ocean[(x_ocean>=start)*(x_ocean<=start+bet)])
    # x_ocean = x_ocean[z_ocean<kk]
    # z_ocean = z_ocean[z_ocean<kk]
    # interplate to the grid length "bet"
    for yy,xx in enumerate(x_grid):
        if len(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])==0:
            continue
        
        if ox[yy-1]< 730:
            oz[yy] = np.min(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
            
        else:
            oz[yy] = np.max(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
        px = xx

    oxx=ox[ox>start]
    oz=oz[ox>start]
    ox=oxx
    return ox,oz

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
    trx = trench_x[frame]
    zslab = w5[xslab>trx]
    xslab = xslab[xslab>trx]
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

            for ele_index in range(len(ele_x[x_ind,:][submantle])):
                maxmax = 9999
                x1 = ele_x[x_ind,:][submantle][ele_index]
                z1 = ele_z[x_ind,:][submantle][ele_index]              
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
    
    Ptotal = np.zeros(len(midslabx))
    dl = np.zeros(len(midslabx))
    torque_length = np.zeros(len(midslabx))
    P0 = np.array((trench_x[frame], trench_z[frame])) * 1000
    Fsux2 = np.zeros(len(xslab))
    Fsux2s = np.zeros(len(xslab))
    beta = np.zeros(len(xslab))
    for ss in range(len(midslabx)-1):
        psub = 0
        ptop = 0
        costheta= 0
        if ss >= 1:
            P3 = np.array((midslabx[ss-1],midslabz[ss-1]))*1e3
            P2 = np.array((xslab[ss],zslab[ss]))*1e3
            P1 = np.array((xslab[ss-1],zslab[ss-1]))*1e3
            dl[ss] = np.linalg.norm(P2-P1)
            torque_length[ss]= np.linalg.norm(P3-P0)
            PPP = np.array([P3-P0,P3-P1])
            costheta = np.dot(PPP[0],PPP[1])/(0.5*dl[ss])/torque_length[ss]
            # print(dl[ss]/1000,torque_length[ss]/1000,PPP/1000,costheta)
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
        beta[ss] = np.arccos(costheta)*180/np.pi
        #beta[ss] = costheta
        Fsu += Ptotal[ss]*dl[ss]*torque_length[ss]*costheta
        Fsux2[ss] = Fsu
        Fsux2s[ss] = Ptotal[ss]*dl[ss]*torque_length[ss]*costheta
        # Fz = (Ptotal*length).sum() # N/m
        # print('Fz=',Fz/1e12)
        # print('Fsu=',Fsu/1e12)
    if frame < 10:
        qq = '00'+str(frame)
    elif frame < 100 and frame >=10:
        qq = '0'+str(frame)
    else:
        qq=str(frame)
    fs.save_4txt(model+'_suctionx_'+str(qq),savepath,xslab[Fsux2!=0],Fsux2[Fsux2!=0],Fsux2s[Fsux2!=0],beta[Fsux2!=0])
    print(savepath+model+'_suctionx_'+str(qq))
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
    # for i in range(55,59):
        loop_time = time.time()
        # fsb[i] = slab_sinking_force(i)
        fsb[i] = slab_sinking_torque(i)
        fsu[i] = suction_force3_torque(i)
        # print(savepath+model+'_suctionx_'+str(i))

    fs.save_3txt(model+'_forces',savepath,fl.time,fsb,fsu)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    time,fbb,fsu = np.loadtxt(savepath+model+'_forces.txt').T
    
    fbb = fd.moving_window_smooth(fbb,8)
    fsu = fd.moving_window_smooth(fsu,8)
    ax.plot(fl.time,fbb,c='#c06c84',label='slab pull (N)',lw=4)
    ax.plot(fl.time,fsu,c="#355c7d",label='suction force (N)',lw=4)
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
