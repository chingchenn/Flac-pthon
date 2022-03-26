#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 13:11:24 2021
@author: jiching
"""
import math
import flac
import os,sys
import numpy as np
import pandas as pd
import gravity as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as f2
import matplotlib.pyplot as plt

#---------------------------------- DO WHAT -----------------------------------
## creat data
vtp                     = 0
trench_location         = 1
dip                     = 1
magma                   = 0
gravity                 = 0
gravity_frame           = 0
melting                 = 0
stack_topo              = 0 
# stack_gem


# plot data
trench_plot             = 1
dip_plot                = 1
magma_plot              = 0
marker_number           = 0
gravity_plot            = 0
phase_plot              = 0
phase_accre             = 0
melting_plot            = 0
force_plot_LR           = 0
force_plot_RF           = 0
vel_plot                = 0
stack_topo_plot         = 0
# stack_gem_plot

#---------------------------------- SETTING -----------------------------------
plt.rcParams["font.family"] = "Times New Roman"
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
sys.path.append("/home/jiching/geoflac/util")
#model = 'k0425'
model = sys.argv[1]
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time=fl.time
#------------------------------------------------------------------------------
def trench(start_vts=1,model_steps=end):
    trench_x=np.zeros(end)
    trench_z=np.zeros(end)
    trench_index=np.zeros(end)
    for i in range(start_vts,model_steps):
        x,z = fl.read_mesh(i)
        sx,sz=f2.get_topo(x,z)
        arc_ind,trench_ind=f2.find_trench_index(z)
        trench_index[i]=trench_ind
        trench_x[i]=sx[trench_ind]
        trench_z[i]=sz[trench_ind]
    return trench_index,trench_x,trench_z

def get_topo(start=1,end_frame=end): 
    topo = [];dis = [];time = []  # do not change to array since the topo database is 3D
    trench_index, xtrench, ztrench = trench(start,end_frame)
    for step in range(start,end_frame):
        x,z = fl.read_mesh(step)
        sx,sz=f2.get_topo(x,z)
        topo.append(sz) 
        dis.append(sx)
        for ii in range(len(sx)):
            time.append(fl.time[step])
    return  dis, time, topo
trenchfile=path+'data/trench_for_'+model+'.csv'
def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
def dip_setting(depth1,depth2):
    depth1=depth1
    depth2=depth2
    return depth1,depth2
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
            crust_x[j] = np.average(ele_x[j,ind_oceanic])
            crust_z[j] = np.average(ele_z[j,ind_oceanic])        
    return crust_x,crust_z
def plate_dip(depth1,depth2):
    angle = np.zeros(end)
    for i in range(1,end):
        crust_x,crust_z = oceanic_slab(i)
        ind_within_80km = (crust_z >= int(depth2)) * (crust_z < int(depth1))
        if not True in (crust_z < int(depth2)):
            continue
        crust_xmin = np.amin(crust_x[ind_within_80km])
        crust_xmax = np.amax(crust_x[ind_within_80km])
        crust_zmin = np.amin(crust_z[ind_within_80km])
        crust_zmax = np.amax(crust_z[ind_within_80km])
        dx = crust_xmax - crust_xmin
        dz = crust_zmax - crust_zmin
        angle[i] = math.degrees(math.atan(dz/dx))
    return time,angle

def plot_phase_in_depth(depth=0):
    time=[];ph=[];xx=[]
    for step in range(end):
        x, z = fl.read_mesh(step+1)
        phase=fl.read_phase(step+1)
        ele_x, ele_z = nodes_to_elements(x,z)
        xt = ele_x[:,0]
        zt = ele_z[:,0]
        pp = np.zeros(xt.shape)
        t = np.zeros(zt.shape)
        t[:]=fl.time[step]
        for gg in range(len(ele_z)):
            pp[gg]=phase[gg,depth]
        time.append(t)
        ph.append(pp)
        xx.append(xt)
    return xx, time, ph
def get_gravity(start=1,end_frame=end):
    fa=[];bg=[];dis=[];time=[];to=[];tom=[]
    for step in range(start+1,end_frame+1):
        px, topo, topomod, fa_gravity, gb_gravity=fg.compute_gravity2(step)
        px *= 10**-3
        topo *= 10**-3
        topomod *=10**-3
        fa_gravity *= 10**5
        gb_gravity *= 10**5
        if gravity_frame:
            if not os.path.exists(savepath+model):
                os.makedirs(savepath+model)
            fs.save_5array('topo-grav.'+str(step), savepath+model, px, topo, topomod, fa_gravity,
                           gb_gravity, 'disX', 'topo', 'topomod', 'free-air', 'bourger')
        fa.append(fa_gravity)
        bg.append(gb_gravity)
        dis.append(px)
        to.append(topo)
        tom.append(topomod)
        for yy in range(len(px)):
            time.append(fl.time[step])
    return dis,time,to,tom,fa,bg
def get_magma(start_vts=1,model_steps=end-1):
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    for i in range(1,end):
        x,z=fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        chamber=fl.read_fmagma(i)
        melt[i] = np.max(mm)
        magma[i] = np.max(chamber)
        arc_vol[i]=np.sum(fl.read_area(i)[phase ==14])/1e6
        yymelt[i]=(fl.read_fmelt(i)*fl.read_area(i)/1e6).sum()
        yychamber[i]=(fl.read_fmagma(i)*fl.read_area(i)/1e6).sum()
    return melt,magma,yymelt,yychamber,arc_vol
magmafile=path+'data/magma_for_'+model+'.csv' 
def count_marker(phase,start=1,end_frame=end):
    mr = np.zeros(end_frame-start)
    for i in range(start,end_frame):
        x,y,age,ph,id=fl.read_markers(i)
        ppp=(ph==phase)
        select1 = id[ppp]
        if len(select1)==0:
            continue
        xk,yk,dk,phk,idk=fl.read_markers(i+1)
        count = 0
        ind_p = idk[(phk==phase)]
        for j in select1:
            if j in ind_p:
                count += 1
        mr[end_frame-i-1]=count   
    return mr
def melting_phase():
    melt_num = np.zeros(end)
    phase_p4=np.zeros(end)
    phase_p9=np.zeros(end)
    phase_p10=np.zeros(end)
    po=np.zeros(end)
    for i in range(1,end):
        c=0;p9=0;p4=0;p10=0
        x, z = fl.read_mesh(i)
        mm=fl.read_fmelt(i)
        phase=fl.read_phase(i)
        for xx in range(len(mm)):
            for zz in range(len(mm[0])):
                if mm[xx,zz] != 0:
                    print()
                    if phase[xx,zz]==9:
                        p9 += 1
                    elif phase[xx,zz]==4:
                        p4 += 1
                    elif phase[xx,zz]==10:
                        p10 += 1
                    c +=1
        pk=c-p4-p9-p10
        melt_num[i]=c
        phase_p4[i]=p4
        phase_p9[i]=p9
        phase_p10[i]=p10
        po[i]=pk
    return phase_p4,phase_p9,phase_p10,po
def get_stack_topo(width=600,ictime=20):
    topo1 = 0;xmean = 0
    fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(xt.shape)
        t[:] = i*0.2
        x_trench = xt[np.argmin(zt)]
        within_plot = (xt>x_trench-width) * (xt<x_trench+width)
        if i >= end-ictime:
            topo1 += zt
            xmean += (xt-x_trench)
        if stack_topo_plot:
            rainbow = cm.get_cmap('gray_r',end)    
            newcolors = rainbow(np.linspace(0, 1, end))
            ax2.plot(xt[within_plot]-x_trench,zt[within_plot],c=newcolors[i])
        xx=(xmean[within_plot]/ictime)
        zz=(topo1[within_plot]/ictime)
    return xx,zz
#------------------------------------------------------------------------------
if vtp:
   file=path+model
   cmd = '''
cd %(file)s 
python /home/jiching/geoflac/util/flacmarker2vtk.py . -1
''' % locals()
   os.system(cmd)
if trench_location:
    print('-----creat trench database-----')
    name='trench_for_'+model
    trench_index,trench_x,trench_z=trench()
    fs.save_3array(name,savepath,time,trench_x,trench_z,
                'time','trench_x','trench_z')
    print('=========== DONE =============')
if dip:
    print('-----creat angle database-----')
    name='plate_dip_of_'+model
    depth1,depth2 = dip_setting(-5,-120)
    time,dip = plate_dip(depth1,depth2)
    fs.save_2array(name,savepath,time,dip,'time','angle')
    print("============ DONE ============")
if gravity:
    print('-----creat gravity database----- ')
    name='gravity_all_'+model
    dis,time,to,tom,fa,bg=get_gravity()
    fs.save_6array(name, savepath, time, dis, to, tom, fa, bg,
                   'time', 'disX', 'topo', 'topomod', 'free-air', 'bourger')
    print('=========== DONE =============')
# if not os.path.exists(magmafile):
if magma:
    print('-----creat magma database-----')
    name='magma_for_'+model
    melt,chamber,yymelt,yychamber,rrr=get_magma()      
    fs.save_6array(name,savepath,time,melt,chamber,yymelt,yychamber,rrr,
                   'time','fmelt','chamber','production','volume','ratio')
    print('=========== DONE =============')
if melting:
    print('-----creat magma database-----')
    name='melting_'+model
    phase_p4,phase_p9,phase_p10,po=melting_phase()
    fs.save_5array(name,savepath,time,phase_p4,phase_p9,phase_p10,po,
                'time','phase_4','phase_9','phase_10','others')
    print('=========== DONE =============')
if stack_topo:
    print('-----creat topo database-----')
    name=model+'_stack_topography'
    xx,zz=get_stack_topo()
    fs.save_2txt(name,savepath,xx,zz)
    print('=========== DONE =============')
##------------------------------------ plot -----------------------------------
if trench_plot:
    print('--- start plotting the trench and topography with time ---')
    name='trench_for_'+model
    #df = pd.read_csv(path+'data/'+name+'.csv')
    df = pd.read_csv(savepath+name+'.csv')
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    dis,time,topo=get_topo(start=1,end_frame=end)
    qqq=ax.scatter(dis,time,c=topo,cmap='gist_earth',vmax=6,vmin=-10)
    cbar=fig.colorbar(qqq,ax=ax)
    ax.plot(df.trench_x[df.trench_x>0],df.time[df.trench_x>0],c='k',lw=2)
    ax.set_xlim(0,dis[-1][-1])
    ax.set_ylim(0,fl.time[-1])
    ax.set_title(str(model)+" Bathymetry Evolution",fontsize=24)
    ax.set_ylabel('Time (Myr)',fontsize=20)
    ax.set_xlabel('Distance (km)',fontsize=20)
    cbar.set_label('Topography (km)',fontsize=20)
    fig.savefig(figpath+model+'_topo.png')
    print('=========== DONE =============')
if dip_plot:
    name = 'plate_dip_of_'+model
    depth1,depth2 = dip_setting(-5,-120)
    df = pd.read_csv(savepath+name+'.csv')
    fig, (ax2)= plt.subplots(1,1,figsize=(10,7))
    ax2.plot(fl.time[df.angle>0],df.angle[df.angle>0],c='royalblue',lw=2)
    ax2.set_xlim(0,fl.time[-1])
    ax2.set_title('Angle Variation of '+str(model),fontsize=24)
    ax2.set_xlabel('Time (Myr)',fontsize=20)
    ax2.set_ylabel('Angel ($^\circ$) from '+str(-depth1)+' to '+str(-depth2)+' depth',fontsize=20)
    ax2.grid()
    fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip.jpg')
if magma_plot:
    name='magma_for_'+model
    df = pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax,ax2,ax3,ax4) = plt.subplots(4,1,figsize=(15,15))
    ax.plot(df.time,df.production,color='tomato')
    ax2.plot(df.time,df.volume,color='orange')
    ax3.bar(df.time,df.fmelt,width=0.1,color='tomato',label='fmelt')
    ax4.bar(df.time,df.chamber,width=0.1,color='orange',label='magma')
    #ax.set_xlabel('Time (Myr)',fontsize=20)
    #ax2.set_xlabel('Time (Myr)',fontsize=20)
    #ax3.set_xlabel('Time (Myr)',fontsize=20)
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('melt * area',fontsize=20)
    ax2.set_ylabel('chamber *area',fontsize=20)
    ax3.set_ylabel('max melt',fontsize=20)
    ax4.set_ylabel('max magma fraction',fontsize=20)
    #ax.set_ylim(0,0.8)
    #ax2.set_ylim(0,10*1e-3)
    #ax3.set_ylim(0,10*1e-3)
    #ax4.set_ylim(0,3*1e-5)
    ax.set_xlim(0,24)
    ax2.set_xlim(0,24)
    ax3.set_xlim(0,24)
    ax4.set_xlim(0,24)
    ax.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax2.tick_params(axis='x', labelsize=16 )
    ax3.tick_params(axis='x', labelsize=16 )
    ax4.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax2.tick_params(axis='y', labelsize=16 )
    ax3.tick_params(axis='y', labelsize=16 )
    ax4.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    fig.savefig(figpath+model+'_magma.png')
#--------------------------------------------------------------------
'''
    fig2,(ax,ax2)=plt.subplots(1,2,figsize=(25,8))
    cb_plot=ax.scatter(df.fmelt,df.chamber,c=df.time,cmap='rainbow')
    ax_cbin =fig2.add_axes([0.13,0.78,0.23,0.03]) 
    cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    ax_cbin.set_title('Myr')
    rrr1=f2.moving_window_smooth(df.ratio,10)
    ax2.plot(df.time,rrr1,color='k',lw=3)
    ax2.plot(df.time,df.ratio,color='gray',linestyle=':')
    ax.set_ylim(0,max(df.chamber))
    ax.set_xlim(0,max(df.fmelt))
    ax.set_ylabel('max magma fraction')
    ax.set_xlabel('max melt fraction')
    ax2.set_xlabel('Myr')
    ax2.set_ylim(0,max(df.ratio))
    ax2.set_xlim(0,max(df.time))
    fig2.savefig(figpath+model+'_ratio.png')
'''
#--------------------------------------------------------------------
if marker_number != 0:
    mr = count_marker(marker_number)
     #plt.plot(mr,c='b')
if gravity_plot:
    name='gravity_for_'+model
    fig, (ax,ax2)= plt.subplots(1,2,figsize=(22,12)) 
    dis,time,to,tom,fa,bg=get_gravity(1,end)
    qqq=ax.scatter(dis,time,c=fa,cmap='Spectral',vmax=400,vmin=-400)
    ax2.scatter(dis,time,c=bg,cmap='Spectral',vmax=400,vmin=-400)
    fig.colorbar(qqq,ax=ax)
    ax2.set_title('bourger gravoty anomaly')
    ax.set_title('free-air gravity anomaly')
    plt.savefig(figpath+model+'_gravity.png')
if phase_plot:
    name = 'phase_for'+model
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    colors = ["#CECCD0","#FF00FF","#8BFF8B","#7158FF","#FF966F",
          "#9F0042","#660000","#524B52","#D14309","#5AB245",
          "#004B00","#008B00","#455E45","#B89FCE","#C97BEA",
          "#525252","#FF0000","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    xt,t,pp= plot_phase_in_depth(depth=0)    
    mmm=ax.scatter(xt,t,c=pp,cmap=phase15,vmin=1, vmax=18)
    ax.set_ylabel("Time (Ma)")
    ax.set_xlabel("Distance (km)")
    ax.set_title(str(model)+" Phase")
    ax.set_ylim(0,t[-1][-1])
    ax.set_xlim(xt[0][0],xt[-1][-1])
    cb_plot1 = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=phase15,vmin=1, vmax=18)
    ax_cbin = fig.add_axes([0.27, 0.03, 0.23, 0.03])
    cb = fig.colorbar(cb_plot1,cax=ax_cbin,orientation='horizontal')
    ax_cbin.set_title('Phase')
    fig.savefig(figpath+model+'_phase.png')
if phase_accre:
    name='trench_for_'+model
    df = pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    dis,time,topo=get_topo(start=1,end_frame=end)
    colors = ["#CECCD0","#FF00FF","#8BFF8B","#7158FF","#FF966F",
          "#9F0042","#660000","#524B52","#D14309","#5AB245",
          "#004B00","#008B00","#455E45","#B89FCE","#C97BEA",
          "#525252","#FF0000","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    xt,t,pp= plot_phase_in_depth(depth=0)
    mmm=ax.scatter(xt,t,c=pp,cmap=phase15,vmin=1, vmax=18)
    ax.set_ylabel("Time (Ma)")
    ax.set_xlabel("Distance (km)")
    ax.set_title(str(model)+" Phase")
    ax.set_ylim(0,t[-1][-1])
    ax.set_xlim(xt[0][0],xt[-1][-1])
    cb_plot1 = ax.scatter([-1],[-1],s=0.1,c=[1],cmap=phase15,vmin=1, vmax=18)
    ax_cbin = fig.add_axes([0.27, 0.03, 0.23, 0.03])
    cb = fig.colorbar(cb_plot1,cax=ax_cbin,orientation='horizontal')
    ax_cbin.set_title('Phase')
    ax.plot(df.trench_x,df.time,c='k',lw=2)    
    fig.savefig(figpath+model+'_acc.png')
if melting_plot:
    name='melting_'+model
    df=pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax) = plt.subplots(1,1,figsize=(18,12))
    ax.bar(df.time,df.phase_p9,width=0.17,color='orange',label='serpentinite ')
    ax.bar(df.time,df.phase_p4,bottom=df.phase_p9,width=0.17,color='seagreen',label='olivine')
    ax.bar(df.time,df.phase_p10,bottom=df.phase_p9+df.phase_p4,width=0.17,color='tomato',label='sediments')
    ax.bar(df.time,df.others,bottom=df.phase_p9+df.phase_p4+df.phase_p10,width=0.17,color='k',label='others')
    ax.set_xlim(0,24)
    ax.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('number of elements in phases',fontsize=20)
    ax.legend(fontsize=25)
    fig.savefig(figpath+model+'_bar_plot_melting.png')
if force_plot_LR:
    filepath = '/home/jiching/geoflac/'+model+'/forc.0'
    fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax.scatter(time,forc_l,c="#4682B4",s=4)
    ax2.scatter(time,forc_r,c="#D2691E",s=4)
    ax.set_xlim(0,time[-1])
    ax.set_title('oceanic side force',fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid()
    ax2.set_xlim(0,time[-1])
    ax2.set_title('continental side force',fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    fig.savefig(figpath+model+'_forc.png')
if force_plot_RF:
    filepath = '/home/jiching/geoflac/'+model+'/forc.0'
    fig2, (ax3)= plt.subplots(1,1,figsize=(10,8))   
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax3.scatter(time,ringforce,c="#483D8B",s=4)
    ax3.set_xlim(0,time[-1])
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax3.grid()
    fig2.savefig(figpath+model+'_ringforc.png')
if vel_plot:
    filepath = '/home/jiching/geoflac/'+model+'/forc.0'
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    fig3, (ax4)= plt.subplots(1,1,figsize=(10,8))
    ax4.plot(time,vl*31545741325,c="#000080",lw=2)
    ax4.set_xlim(0,time[-1])
    ax4.set_title('oceanic side velocity',fontsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.grid()
    ax4.set_xlabel('Time (Myr)',fontsize=16)
    ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
    fig3.savefig(figpath+model+'_vel.png')

if stack_topo_plot:
    name=model+'_stack_topography.txt'
    xx,zz=get_stack_topo()
    xmean,ztop=np.loadtxt(path+'data/'+name).T
    ax2.plot(xx,zz,c="#000080",lw=3)
    fig2.savefig(figpath+model+'_topo_analysis.png')