#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:42:12 2022

@author: ji-chingchen
"""

import sys, os
import numpy as np
import flac
import matplotlib
from matplotlib import cm
import function_for_flac as fd
import function_savedata as fs
from scipy import interpolate
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.figsize"] = (10,12)
model = sys.argv[1]
#model = 'b0702k'
#frame = int(sys.argv[2])
path='/home/jiching/geoflac/'
#path='/Users/ji-chingchen/Desktop/model/'
path = '/scratch2/jiching/22summer/'
#path = 'D:/model/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Users/ji-chingchen/Desktop/data/'
#savepath = 'D:/model/data/'
figpath='/home/jiching/geoflac/figure/'
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
    Fsb = 0; g = 10
    for z_ind in range(1,len(ele_z[0])):
    #for z_ind in range(11,15):
        ind_eclogite = (phase[:,z_ind]==phase_eclogite) + (phase[:,z_ind] == phase_eclogite_1) + (phase[:,z_ind] == phase_hydratedmantle) + (phase[:,z_ind] == phase_mantle2)
        man_eclogite = (phase[:,z_ind]== phase_mantle1) + (phase[:,z_ind] == phase_serpentinite)
        if True in ind_eclogite and True in man_eclogite:
            den_mantle = np.average(density[man_eclogite,z_ind])
            den_eco = np.average(density[ind_eclogite,z_ind])
            rho_diff[z_ind] = den_eco - den_mantle
            slab_area[z_ind] = area[ind_eclogite,z_ind].sum()
        Fsb += rho_diff[z_ind] * g * slab_area[z_ind]
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
###---------------------- Mantle suction force with time -------------------------------
#for frame in range(1,end):
frame = 54
x, z = fl.read_mesh(frame)
ele_x, ele_z = flac.elem_coord(x, z)
phase = fl.read_phase(frame)
density = fl.read_density(frame)
pressure = fl.read_pres(frame)
area = fl.read_area(frame)

aatop = 0;aasub = 0
Fsub = 0; Ftop = 0
onepre=pressure.flatten()
a,b=np.polyfit(onepre,ele_z.flatten(),deg=1)
fit=(ele_z.flatten()-b)/a
dypre=(onepre-fit).reshape(len(phase),len(phase[0]))*1e8  # N/m^2
#static = -density.flatten() * 10*1000*ele_z.flatten() # N/m^2
#dypre = (pressure.flatten()*1e8+static).reshape(len(phase),len(phase[0]))/1e6 #MPa

# colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
#           "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
#           "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
#           "#F67280","#00FF00","#FFFF00","#7158FF"]
# phase15= matplotlib.colors.ListedColormap(colors)
# plt.scatter(ele_x,ele_z,c = phase,cmap = phase15,vmin = 1,vmax = 20, s= 40)
# plt.ylim(-450,-0)
# plt.xlim(500,800)
# plt.axes().set_aspect('equal')
final_ind = int(trench_index[frame])
ind_trench = int(trench_index[frame])
xoceanic = np.zeros(len(ele_z)-ind_trench)
#ffsub = np.zeros(len(ele_z)-ind_trench)
#fftop = np.zeros(len(ele_z)-ind_trench)
#asub = np.zeros(len(ele_z)-ind_trench)
#atop = np.zeros(len(ele_z)-ind_trench)
iiind = np.zeros(len(ele_z)-ind_trench)
length = 0
for ii,x_ind in enumerate(range(ind_trench,len(ele_z))):
#for ii,x_ind in enumerate(range(196,210)):
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
            #asub[ii] = (area[x_ind,av_oc_ind:][submantle]).sum()
            #atop[ii] = (area[x_ind,:av_oc_ind][topmantle]).sum()
            
            #ffsub[ii] = ((dypre*area)[x_ind,av_oc_ind:][submantle]).sum()
            #fftop[ii] = ((dypre*area)[x_ind,:av_oc_ind][topmantle]).sum()
            Fsub += ((dypre*area)[x_ind,av_oc_ind:][submantle]).sum()
            Ftop += ((dypre*area)[x_ind,:av_oc_ind][topmantle]).sum()
            # plt.scatter(ele_x[x_ind,av_oc_ind:][submantle],\
            #             ele_z[x_ind,av_oc_ind:][submantle],\
            #             c= 'b',s = 5)
            # plt.scatter(ele_x[x_ind,:av_oc_ind][topmantle],\
            #             ele_z[x_ind,:av_oc_ind][topmantle],\
            #             c= 'r',s = 5)  
            # plt.scatter(ele_x[x_ind,av_oc_ind],ele_z[x_ind,av_oc_ind],c = 'k') 
            
            if ii > 1:
                dx = ele_x[x_ind,av_oc_ind]-ele_x[x_ind-1,int(iiind[ii-1])]
                dz = ele_z[x_ind,av_oc_ind]-ele_z[x_ind-1,int(iiind[ii-1])]
                length += np.sqrt(dx**2 + dz**2)*1000
if length != 0 and aasub !=0 and aatop !=0:
    Fsu = (Fsub/aasub-Ftop/aatop)*length
else:
    Fsu = 0

print(Fsu)  
   
    
def suction_force(frame):
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = flac.elem_coord(x, z)
    phase = fl.read_phase(frame)
    pressure = fl.read_pres(frame)
    area = fl.read_area(frame)
    
    aatop = 0;aasub = 0
    Fsub = 0; Ftop = 0
    onepre=pressure.flatten()
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
    fsb=np.zeros(end)
    ft=np.zeros(end)
    fsu=np.zeros(end)
    ratio=np.zeros(end)
    for i in range(1,end):
        fsb[i] = slab_sinking_force(i)
        ft[i] = mantle_traction_force(i)
        fsu[i] = suction_force(i)
        if ft[i] ==0:
            ratio[i] = 0
        else:
            ratio[i] = fsb[i]/ft[i]

    fs.save_5txt(model+'_forces','/home/jiching/geoflac/data/',fl.time,fsb,ft,fsu,ratio)
    fig, (ax)= plt.subplots(1,1,figsize=(10,6))
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
    ax.set_yscale('log')
    ax.set_title('Forces of '+model,fontsize=20)
    fig.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force.png')
    # fig2, (ax2)= plt.subplots(1,1,figsize=(10,6))
    # ratio_f = fd.moving_window_smooth(ratio[ratio>0],5)
    # ax2.plot(fl.time[ratio>0],ratio_f,c="#355c7d",label='ratio of these forces)',lw=4)
    # ax2.set_xlabel('Time (Myr)',fontsize=16)
    # ax2.set_xlim(0, fl.time[-1])
    # ax2.tick_params(axis='x', labelsize=16)
    # ax2.tick_params(axis='y', labelsize=16)
    # ax2.grid()
    # ax2.spines['bottom'].set_linewidth(bwith)
    # ax2.spines['top'].set_linewidth(bwith)
    # ax2.spines['right'].set_linewidth(bwith)
    # ax2.spines['left'].set_linewidth(bwith)
    #fig2.savefig('/home/jiching/geoflac/figure/'+model+'_slab_force_ratio.png')
