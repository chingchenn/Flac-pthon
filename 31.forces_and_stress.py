#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

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
model = 'b0804m'
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
###----------------------- Slab sinking force with time-------------------------------
#    Fsb = (rho_mantle-rho_slab)(z) * g * area_of_slab 
def slab_sinking_force(frame):
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    density = fl.read_density(frame)
    area = fl.read_area(frame)
    rho_diff = np.zeros(nez)
    slab_area = np.zeros(nez)
    temp = fl.read_temperature(frame)
    ele_temp = temp_elements(temp)
    Fsb = 0; g = 10
    for z_ind in range(1,len(ele_z[0])):
        ind_eclogite = (phase[:,z_ind]==phase_eclogite) + (phase[:,z_ind] == phase_eclogite_1) #+ (phase[:,z_ind] == phase_hydratedmantle) + (phase[:,z_ind] == phase_mantle2)
        man_eclogite = (phase[:,z_ind]== phase_mantle1) + (phase[:,z_ind] == phase_serpentinite)+ (phase[:,z_ind] == phase_hydratedmantle)#+ (phase[:,z_ind] == phase_schist)
        ind_slab = (ele_temp[:,z_ind]<700)*((phase[:,z_ind]== phase_mantle1)+ (phase[:,z_ind]==phase_eclogite)+ (phase[:,z_ind] == phase_eclogite_1)+ (phase[:,z_ind] == phase_hydratedmantle))
        #print(man_eclogite)
        if True in ind_eclogite and True in man_eclogite:
            den_mantle = np.average(density[man_eclogite,z_ind])
            den_eco = np.average(density[ind_eclogite,z_ind])
            if den_eco < den_mantle:
                #print(frame,z_ind, den_eco-den_mantle)
                continue
            rho_diff[z_ind] = den_eco - den_mantle
            slab_area[z_ind] = area[ind_slab,z_ind].sum()
            #if True in ind_slab:
                #print(ele_x[ind_slab],z_ind)
        Fsb += rho_diff[z_ind] * g * slab_area[z_ind]
    print('Fsb=',Fsb/1e12)
    return Fsb # N/m (2D)

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
            #print(qq,last_deep)
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

# def find_slab_median_index(frame):
#     final_ind = int(trench_index[frame])
#     ind_trench = int(trench_index[frame])
#     xoceanic = np.zeros(len(ele_z)-ind_trench)
#     slab_x = np.zeros(len(ele_z)-ind_trench)
#     slab_z = np.zeros(len(ele_z)-ind_trench)
#     for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
#     # for ii,x_ind in enumerate(range(140,272)):
#         ind_oceanic = (phase[x_ind,:] == phase_oceanic) + (phase[x_ind,:] == phase_eclogite)+(phase[x_ind,:] == phase_oceanic_1) + (phase[x_ind,:] == phase_eclogite_1)
#         subducted_sed = (phase[x_ind,:] == phase_sediment) + (phase[x_ind,:] ==phase_sediment_1) + (phase[x_ind,:] == phase_schist)
#         if True in (ind_oceanic+subducted_sed):
#             if (x_ind - final_ind) > 2:
#                 break
#             final_ind = x_ind
    
#             oceanic_plate_index = [ww for ww, x in enumerate(ind_oceanic+subducted_sed) if x]
#             xoceanic[ii] = int(np.median(oceanic_plate_index))
#             xstd = np.std(oceanic_plate_index)
#             oo = oceanic_plate_index
#             if xstd > 15:
#                 print(x_ind)
#                 oo = np.array(oceanic_plate_index)[oceanic_plate_index>=xoceanic[ii-1]]
#                 if len(oo)==0:
#                     continue
#                 xoceanic[ii] = int(np.median(oo))
#             av_oc_ind = int(np.median(oo))
#             slab_x[ii] = ele_x[x_ind,av_oc_ind]
#             slab_z[ii] = ele_z[x_ind,av_oc_ind]
#         # print(ii,av_oc_ind,oo,x_ind)
#     return slab_x,slab_z

def dynamic_pressure(frame):
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    pressure = -fl.read_pres(frame)
    onepre=pressure.flatten()
    a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
    fit=(ele_z.flatten()-b)/a
    dypre=(onepre-fit).reshape(len(ele_x),len(ele_x[0]))*1e8  # N/m^2
    return dypre # N/m^2
###---------------------- Mantle suction force with time -------------------------------

# depth1 = -350
# depth2 = -10
# frame = 30
# x, z = fl.read_mesh(frame)
# ele_x, ele_z = flac.elem_coord(x, z)
# phase = fl.read_phase(frame)
# dpre = dynamic_pressure(frame) # N/m^2
# slab_x,slab_z = find_slab_median_index2(frame)

# xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
# zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]


# z1=np.polyfit(xslab,zslab,4)
# p4=np.poly1d(z1)
# w1=p4(xslab) # 4st poly 
# p3=np.polyder(p4,1)
# p2=np.polyder(p4,2)
# w2=p3(xslab) # 3st poly 
# w3=p2(xslab) # 2st poly

# z5 =  np.polyfit(xslab,zslab,5)
# p5 = np.poly1d(z5)
# w5 = p5(xslab) # 5st poly

# list_subslab =  [[] for i in range(len(zslab))]
# list_topslab =  [[] for i in range(len(zslab))]
# ind_trench = int(trench_index[frame])

# colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
#           "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
#           "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
#           "#F67280","#00FF00","#FFFF00","#7158FF"]
# phase15= matplotlib.colors.ListedColormap(colors)
# cm = plt.cm.get_cmap('RdYlBu_r')
# plt.scatter(ele_x,ele_z,c=-dpre/1e6,cmap=cm,vmin=-200, vmax=200,s=40)
# # plt.scatter(ele_x,ele_z,c=phase,cmap=phase15,vmin=1, vmax=20,s=100)
# plt.ylim(-300,-0)
# plt.xlim(300,1000)
# plt.axes().set_aspect('equal')
# #plt.scatter(xslab,zslab,c = 'k',s = 50)
# plt.scatter(xslab,w5,c = 'k',s = 30)
# zslab = w5

# #p = 0
# for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
# #for ii,x_ind in enumerate(range(ind_trench,170)):
#     zind = np.where(ele_z[x_ind,:]==slab_z[ii])[0]
#     if len(zind)==0:
#         zind = 0
#     else:
#         zind = int(zind)
#     # Choose the submantle area
#     isabove = lambda p,a,b : np.cross(p-a,b-a)<0
#     submantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*((phase[x_ind,:]==phase_mantle1) + \
#             (phase[x_ind,:]==phase_mantle2) +(phase[x_ind,:]==phase_hydratedmantle))
#     if True in submantle:
#         for ele_index in range(len(ele_x[x_ind,:][submantle])):
#             maxmax = 9999
#             x1 = ele_x[x_ind,:][submantle][ele_index]
#             z1 = ele_z[x_ind,:][submantle][ele_index]
#             for slab_ind in range(len(xslab)):
#                 dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
#                 if dis<maxmax:
#                     maxmax = dis
#                     choosex = xslab[slab_ind]
#                     choosez = zslab[slab_ind]
#                     chooseind = slab_ind
#             if maxmax < 50:
#                 if chooseind  == len(xslab)-1:
#                     kk = len(xslab)-1
#                 else:
#                     kk = chooseind+1
#                 pointp = np.array([x1,z1])
#                 pointa = np.array([choosex,choosez])
#                 pointb = np.array([xslab[kk],zslab[kk]])
#                 if isabove(pointp,pointa,pointb)==False:
#                     z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
#                     list_subslab[chooseind].append([x_ind,z_ind])
#                     #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='b',s = 5)
#                 #elif isabove(pointp,pointa,pointb)==True:
#                     #z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
#                     #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 
#     #Choose the top mantle area
#     if slab_z[ii]==0:
#         slab_z[ii]=-410
#     #topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*\
#     #((phase[x_ind,:]==phase_mantle1) + (phase[x_ind,:]==phase_mantle2) + \
#       #       (phase[x_ind,:]==phase_serpentinite)+(ele_x[x_ind,:]>max(slab_x)))*(ele_z[x_ind,:]>slab_z[ii])
#     topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*\
#     ((phase[x_ind,:]==phase_mantle1) + (phase[x_ind,:]==phase_mantle2) + \
#             (phase[x_ind,:]==phase_serpentinite))
#     if True in topmantle:
#         for ele_index in range(len(ele_x[x_ind,:][topmantle])):
#             maxmax = 9999
#             x1 = ele_x[x_ind,:][topmantle][ele_index]
#             z1 = ele_z[x_ind,:][topmantle][ele_index]
#             for slab_ind in range(len(xslab)):
#                 dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
#                 if dis<maxmax:
#                     maxmax = dis
#                     choosex = xslab[slab_ind]
#                     choosez = zslab[slab_ind]
#                     chooseind = slab_ind
#             if maxmax < 50:
#                 if chooseind  == len(xslab)-1:
#                     kk = len(xslab)-1
#                 else:
#                     kk = chooseind+1
#                 pointp = np.array([x1,z1])
#                 pointa = np.array([choosex,choosez])
#                 pointb = np.array([xslab[kk],zslab[kk]])
#                 if isabove(pointp,pointa,pointb)==True:
#                     z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
#                     list_topslab[chooseind].append([x_ind,z_ind])
#                     #plt.scatter(ele_x[x_ind,z_ind],ele_z[x_ind,z_ind],c='r',s = 5) 

# Ptotal = np.zeros(len(xslab))
# dl = np.zeros(len(xslab))
# length = np.zeros(len(xslab))
# for ss in range(len(xslab)-1):
# #for ss in range(139,140):
#     psub = 0
#     ptop = 0
#     if ss >= 1:
#         x0 = xslab[ss]
#         z0 = zslab[ss]
#         x1 = xslab[ss-1]
#         z1 = zslab[ss-1]
#         dl[ss] = np.sqrt((-(zslab[ss]-zslab[ss-1]))**2+(-(xslab[ss]-xslab[ss-1]))**2)*1e3  # m
#         vec = [x0-x1,z0-z1]
#         ref = [1,0]
#         #print(vec)
#         if np.linalg.norm(vec)==0:
#             print(ss,x1,z1)
#             continue
#         #sinan = np.cross(vec,ref)/np.linalg.norm(vec)
#         cosan = np.dot(vec,ref)/np.linalg.norm(vec)
#         length[ss] = dl[ss]*cosan
        
#     pres = np.zeros(len(list_subslab[ss]))
#     pret = np.zeros(len(list_topslab[ss]))
#     for rr in range(len(list_subslab[ss])):
#         indx = list_subslab[ss][rr][0]
#         indz = list_subslab[ss][rr][1]
#         xm = ele_x[indx,indz]
#         zm = ele_z[indx,indz]
#         #plt.scatter(xm,zm,c='b',s = 20)
#         pres[rr] = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
#     if len(list_subslab[ss]) !=0:
#         psub = np.average(pres)
#     for rr in range(len(list_topslab[ss])):
#         indx = list_topslab[ss][rr][0]
#         indz = list_topslab[ss][rr][1]
#         xm = ele_x[indx,indz]
#         zm = ele_z[indx,indz]
#         #plt.scatter(xm,zm,c='r',s = 20)
#         pret[rr] = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
#     if len(list_topslab[ss]) !=0:
#         ptop = np.average(pret)
#     Ptotal[ss] = psub-ptop
# Fz = (Ptotal*length).sum()
# print('Fz=',Fz/1e13)

    
# Fsu = (Ptotal*dz).sum() # N/m

def suction_force2(frame):
    depth1 = -350
    depth2 = -40
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    dpre = dynamic_pressure(frame) # N/m^2
    slab_x,slab_z = find_slab_median_index2(frame)
    xslab = slab_x[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    zslab = slab_z[(slab_x>0)*(slab_z>depth1)*(slab_z<depth2)]
    list_subslab =  [[] for i in range(len(zslab))]
    list_topslab =  [[] for i in range(len(zslab))]
    ind_trench = int(trench_index[frame])
    
    for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
        zind = np.where(ele_z[x_ind,:]==slab_z[ii])[0]
        if len(zind)==0:
            zind = 0
        else:
            zind = int(zind)
        # Choose the submantle area
        submantle = (ele_z[x_ind,zind:]< depth2)*(ele_z[x_ind,zind:]> depth1)*((phase[x_ind,zind:]==phase_mantle1) + \
                (phase[x_ind,zind:]==phase_mantle2) + \
                (phase[x_ind,zind:]==phase_serpentinite)+\
                (phase[x_ind,zind:]==phase_hydratedmantle))*(ele_x[x_ind,zind:]<max(slab_x))
        if True in submantle:
            for ele_index in range(len(ele_x[x_ind,zind:][submantle])):
                maxmax = 9999
                x1 = ele_x[x_ind,zind:][submantle][ele_index]
                z1 = ele_z[x_ind,zind:][submantle][ele_index]
                for slab_ind in range(len(xslab)):
                    dis = np.sqrt((xslab[slab_ind]-x1)**2+(zslab[slab_ind]-z1)**2)
                    if dis<maxmax:
                        maxmax = dis
                        choosex = xslab[slab_ind]
                        choosez = zslab[slab_ind]
                        chooseind = slab_ind
                if maxmax < 50:
                    z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                    list_subslab[chooseind].append([x_ind,z_ind])
        # Choose the top mantle area
        if slab_z[ii]==0:
            slab_z[ii]=-410
        topmantle = (ele_z[x_ind,:]< depth2)*(ele_z[x_ind,:]> depth1)*((phase[x_ind,:]==phase_mantle1) + \
                (phase[x_ind,:]==phase_mantle2) + \
                (phase[x_ind,:]==phase_serpentinite)+(ele_x[x_ind,:]>max(slab_x)))*(ele_z[x_ind,:]>slab_z[ii])
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
                    z_ind = np.where(ele_z[x_ind,:]==z1)[0][0]
                    list_topslab[chooseind].append([x_ind,z_ind])
    
    Ptotal = np.zeros(len(xslab))
    dz = np.zeros(len(xslab))
    for ss in range(len(xslab)):
        psub = 0
        ptop = 0
        x0 = xslab[ss]
        z0 = zslab[ss]
        for rr in range(len(list_subslab[ss])):
            indx = list_subslab[ss][rr][0]
            indz = list_subslab[ss][rr][1]
            x1 = ele_x[indx,indz]
            z1 = ele_z[indx,indz]
            vec = [x0-x1,z0-z1]
            ref = [1,0]
            if np.linalg.norm(vec)==0:
                continue
            sinan = np.cross(vec,ref)/np.linalg.norm(vec)
            pre = dpre[list_subslab[ss][rr][0],list_subslab[ss][rr][1]] # N/m^2
            pre_vetical = pre*sinan
            psub += pre_vetical
    
        for rr in range(len(list_topslab[ss])):
            indx = list_topslab[ss][rr][0]
            indz = list_topslab[ss][rr][1]
            x1 = ele_x[indx,indz]
            z1 = ele_z[indx,indz]
            vec = [x0-x1,z0-z1]
            ref = [1,0]
            sinan = np.cross(vec,ref)/np.linalg.norm(vec)
            pre = dpre[list_topslab[ss][rr][0],list_topslab[ss][rr][1]] # N/m^2
            pre_vetical = pre*sinan
            ptop += pre_vetical
        Ptotal[ss] = psub-ptop
    
        if ss >= 1:
            dz[ss] = -(zslab[ss]-zslab[ss-1])*1e3
    Fsu = (Ptotal*dz).sum() # N/m
    return Fsu
    

def suction_force3(frame):
    depth1 = -350
    depth2 = -10
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    dpre = dynamic_pressure(frame) # N/m^2
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
#if __name__ == '__kk__':
    fsb=np.zeros(end)
    ft=np.zeros(end)
    fsu=np.zeros(end)
    ratio=np.zeros(end)
    for i in range(5,end):
        fsb[i] = slab_sinking_force(i)
        ft[i] = mantle_traction_force(i)
        #fsu[i] = suction_force(i)
        fsu[i] = suction_force3(i)
    #    if fsu[i] ==0:
    #        ratio[i] = 0
    #    else:
    #        ratio[i] = fsb[i]/fsu[i]

    fs.save_5txt(model+'_forces',savepath,fl.time,fsb,ft,fsu,ratio)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
    fl.time,fbb,tt,fsu,ratio = np.loadtxt(savepath+model+'_forces.txt').T
    
    sb = fd.moving_window_smooth(fsb[fsb>0],8)
    tt = fd.moving_window_smooth(fsu[fsu>0],8)
    ax.plot(fl.time[fsb>0],sb,c='#c06c84',label='slab pull (N/m)',lw=4)
    ax.plot(fl.time[fsu>0],tt,c="#355c7d",label='suction force (N/m)',lw=4)
    #ax.scatter(fl.time[fsb>0],fsb[fsb>0],c='#c06c84',label='slab pull (N/m)')
    #ax.scatter(fl.time[ft>0],ft[ft>0],c="#355c7d",label='traction force (N/m)')
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
    fig.savefig(figpath+model+'_slab_force.png')
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
