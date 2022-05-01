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
import function_for_flac as fd
import matplotlib.pyplot as plt
from numpy import unravel_index

#---------------------------------- DO WHAT -----------------------------------
## creat data
vtp                     = 0
trench_location         = 0
dip                     = 0
magma                   = 0
melting_loc             = 0
gravity                 = 0
gravity_frame           = 0
melting                 = 0
stack_topo              = 0
stack_gem	        = 1
flat_duraton            = 0

# plot data
trench_plot             = 0
dip_plot                = 0
magma_plot              = 0
metloc_plot  	    	= 0
marker_number           = 0
gravity_plot            = 0
phase_plot              = 0
phase_accre             = 0
melting_plot            = 0
force_plot_LR           = 0
force_plot_RF           = 0
vel_plot                = 0
stack_topo_plot         = 0
stack_gem_plot	    	= 0
flat_slab_plot          = 0

#---------------------------------- SETTING -----------------------------------
path = '/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
savepath='/home/jiching/geoflac/data/'
figpath='/home/jiching/geoflac/figure/'
model = sys.argv[1]
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
time = fl.time
#------------------------------------------------------------------------------
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
def get_topo(start=1,end_frame=end): 
    topo = [];dis = [];time = []  # do not change to array since the topo database is 3D
    trench_index, xtrench, ztrench = trench(start,end_frame)
    for step in range(start,end_frame):
        x,z = fl.read_mesh(step)
        sx,sz=fd.get_topo(x,z)
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
            crust_x[j] = np.max(ele_x[j,ind_oceanic])
            crust_z[j] = np.max(ele_z[j,ind_oceanic])        
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
        if dz ==0:
            continue
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
def melting_location(start_vts=1,model_steps=end-1):
    melt=np.zeros(end)
    x_melt=np.zeros(end)
    for i in range(1,end):
        x,z=fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        melt[i] = np.max(mm)
        maxindex=unravel_index(mm.argmax(),mm.shape)[0]
        ele_x, ele_z = nodes_to_elements(x,z)
        tren_x = (ele_z[:,0]).argmin()
        x_melt[i] = ele_x[maxindex,0]-ele_x[tren_x,0]
    return melt,x_melt,time
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
    phase_p3=np.zeros(end)
    phase_p4=np.zeros(end)
    phase_p9=np.zeros(end)
    phase_p10=np.zeros(end)
    for i in range(1,end):
        c=0;p9=0;p4=0;p10=0;p3=0
        mm=fl.read_fmelt(i)
        phase=fl.read_phase(i)
        area = fl.read_area(i)
        for xx in range(len(mm)):
            for zz in range(len(mm[0])):
                if mm[xx,zz] != 0:
                    if phase[xx,zz]==9:
                        p9 += area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==4:
                        p4 +=area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==10 or phase[xx,zz]==5:
                        p10 += area[xx,zz]*mm[xx,zz]/1e6
                    elif phase[xx,zz]==3:
                        p3 += area[xx,zz]*mm[xx,zz]/1e6
                    c +=1
        pk=c-p4-p9-p10-p3
        melt_num[i]=c
        phase_p3[i]=p3
        phase_p4[i]=p4
        phase_p9[i]=p9
        phase_p10[i]=p10
    return fl.time,phase_p3,phase_p4,phase_p9,phase_p10
def get_stack_topo(width=600,ictime=20):
    topo1 = 0;xmean = 0
    for i in range(end-ictime,end):
        x, z = fl.read_mesh(i)
        xt = x[:,0]
        zt = z[:,0]
        t = np.zeros(xt.shape)
        t[:] = i*0.2
        x_trench = xt[np.argmin(zt)]
        within_plot = (xt>x_trench-width) * (xt<x_trench+width)
        topo1 += zt
        xmean += (xt-x_trench)
        finx = (xt-x_trench)
        finz = zt
    xx=xmean[within_plot]/ictime
    zz=topo1[within_plot]/ictime
    return xx,zz,finx,finz
def get_stack_geometry(ictime=20,width=700):
    stslab = 0;xmean=0
    for i in range(end-ictime,end):
        crust_x,crust_z = oceanic_slab(i)
        x, z = fl.read_mesh(i)
        ele_x, ele_z = nodes_to_elements(x,z)
        x_trench = ele_x[:,0][np.argmin(ele_z[:,0])]
        within_plot = (ele_x[:,0]>x_trench-width)* (crust_z < 0)
        stslab += crust_z
        xmean += (crust_x-x_trench)
        finx = crust_x-x_trench
        finz = crust_z
    xx=xmean[within_plot]/ictime
    zz=stslab[within_plot]/ictime
    return xx[xx>0][:-1],zz[xx>0][:-1],finx,finz
def flat_slab_duration():
    phase_oceanic = 3;phase_ecolgite = 13
    bet = 2;find_flat_dz1=[];find_flat_dz2=[];flat_slab_length=[];flat_slab_depth=[]
    for i in range(1,end):
        x, z = fl.read_mesh(i)
        mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(i)  
        trench_ind = np.argmin(z[:,0]) 
        x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]
        x_ocean = mx[(phase==phase_ecolgite)+(phase==phase_oceanic)]
        z_ocean = mz[(phase==phase_ecolgite)+(phase==phase_oceanic)]
        if z_trench> -2 or min(z_ocean)>-200:
            continue
        start = math.floor(x_trench)
        final = math.floor(np.max(x_ocean))
        x_grid = np.arange(start,final,bet)
        ox = np.zeros(len(x_grid))
        oz = np.zeros(len(x_grid))
        px = start-bet
        kk=np.max(z_ocean[(x_ocean>=start) *(x_ocean<=start+bet)])
        x_ocean = x_ocean[z_ocean<kk]
        z_ocean = z_ocean[z_ocean<kk]
        for yy,xx in enumerate(x_grid):
            if len(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])==0:
                continue
            oz[yy] = np.average(z_ocean[(x_ocean>=px)*(x_ocean<=xx)])
            ox[yy] = np.average(x_ocean[(x_ocean>=px)*(x_ocean<=xx)])
            px = xx
        oxx=ox[ox>start]
        oz=oz[ox>start]
        ox=oxx
    ### =========================== polynomial ===========================
        z1=np.polyfit(ox,oz,4)
        p4=np.poly1d(z1)
        w1=p4(ox)
        p3=np.polyder(p4,1)
        p2=np.polyder(p4,2)
        w2=p3(ox)
        w3=p2(ox)
        cc=-1;ff1=[]
        for rr,oo in enumerate(w2):
            if cc*(oo+0.2)<0:
                ff1.append(ox[rr])
            cc = oo+0.2
        mm=-1;ff2=[]
        for pp,uu in enumerate(w3):
            if mm*uu<0:
                ff2.append(ox[pp])
            mm = uu  
        if len(ff2)>1 and (ff2[1]-ff2[0])>80 and ff2[0]>start:
            find_flat_dz2.append(fl.time[i])
            if len(ff1)>2 and (ff1[-1]-ff1[-2])>50:
                find_flat_dz1.append(fl.time[i])
                flat_slab_length.append(ff1[-1]-ff1[-2])
                depth=np.average(w1[(ox>=ff1[-2])*(ox<ff1[-1])])
                flat_slab_depth.append(depth)
    return find_flat_dz2,find_flat_dz1,flat_slab_length,flat_slab_depth
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
    fs.save_3array(name,savepath,fl.time,trench_x,trench_z,
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
if magma:
    print('-----creat magma database-----')
    name='magma_for_'+model
    melt,chamber,yymelt,yychamber,rrr=get_magma()
    fs.save_5txt(name,savepath,melt,chamber,yymelt,yychamber,rrr)
    print('=========== DONE =============')
if melting_loc:
    print('-----creat melt database-----')
    name='metloc_for_'+model
    melt,xmelt,time=melting_location()
    fs.save_3txt(name,savepath,time,melt,xmelt)
    print('=========== DONE =============')
if melting:
    print('-----creat melting database-----')
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p9,phase_p10=melting_phase()
#    fs.save_5array(name,savepath,time,phase_p4,phase_p9,phase_p10,po,
#                'time','phase_4','phase_9','phase_10','others')
    fs.save_5txt(name,savepath,time,phase_p3,phase_p4,phase_p9,phase_p10)
    print('=========== DONE =============')
if stack_topo:
    print('-----creat topo database-----')
    name=model+'_stack_topography'
    xx,zz,fx,fz=get_stack_topo()
    fs.save_2txt(name,savepath,xx,zz)
    fs.save_2txt(model+'_final_topography',savepath,fx,fz)
    print('=========== DONE =============')
if stack_gem:
    print('-----creat geometry database-----')
    name=model+'_stack_slab'
    xx,zz,fx,fz=get_stack_geometry()
    fs.save_2txt(name,savepath,xx,zz)
    fs.save_2txt(model+'_final_slab',savepath,fx,fz)
    print('=========== DONE =============')
if flat_duraton:
    print('-----creat flattime database-----')
    name=model+'_flatslab_duration2'
    f2,f1,length,dep=flat_slab_duration()
    fs.save_1txt(name,savepath,f2)
    name=model+'_flatslab_time_len'
    fs.save_3txt(name,savepath,f1,length,dep)
    print('=========== DONE =============')
##------------------------------------ plot -----------------------------------
if trench_plot:
    print('----- plotting topography-----')
    name='trench_for_'+model
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
    print('------plotting dip angle------')
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
    bwith=3
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip.jpg')
    print('=========== DONE =============')
if magma_plot:
    print('--------plotting magma--------')
    name = 'magma_for_'+model+'.txt'
    filepath = '/home/jiching/geoflac/data/'
    temp1 = np.loadtxt(filepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    fig, (ax,ax2,ax3,ax4) = plt.subplots(4,1,figsize=(15,15))
    ax.plot(fl.time,yymelt,color='tomato')
    ax2.plot(fl.time,yychamber,color='orange')
    ax3.bar(fl.time,melt,width=0.1,color='tomato',label='fmelt')
    ax4.bar(fl.time,chamber,width=0.1,color='orange',label='magma')
    ax4.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('melt * area',fontsize=20)
    ax2.set_ylabel('chamber *area',fontsize=20)
    ax3.set_ylabel('max melt',fontsize=20)
    ax4.set_ylabel('max magma fraction',fontsize=20)
    ax.set_xlim(0,fl.time[-1]);ax.grid()
    ax2.set_xlim(0,fl.time[-1]);ax2.grid()
    ax3.set_xlim(0,fl.time[-1]);ax3.grid()
    ax4.set_xlim(0,fl.time[-1]);ax4.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax2.tick_params(axis='x', labelsize=16 )
    ax3.tick_params(axis='x', labelsize=16 )
    ax4.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax2.tick_params(axis='y', labelsize=16 )
    ax3.tick_params(axis='y', labelsize=16 )
    ax4.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    ax3.spines['bottom'].set_linewidth(bwith)
    ax3.spines['top'].set_linewidth(bwith)
    ax3.spines['right'].set_linewidth(bwith)
    ax3.spines['left'].set_linewidth(bwith)
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_magma.png')
    print('=========== DONE =============')
if metloc_plot:
    print('-------plotting gravity-------')
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    time,melt,xmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    qqq=ax.scatter(time[melt>0],xmelt[melt>0],c=melt[melt>0],cmap='OrRd',vmin=0.0,vmax=0.05)
    cbar=fig.colorbar(qqq,ax=ax)
    ax.set_ylim(0,400)
    ax.set_xlim(0,30)
    ax.set_title(str(model)+" Melting location",fontsize=24)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    cbar.set_label('Melting %',fontsize=20)
    ax.set_ylabel('Distance with trench (km)',fontsize=20)
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_metloc.png')
    print('=========== DONE =============')
#--------------------------------------------------------------------
'''
    fig2,(ax,ax2)=plt.subplots(1,2,figsize=(25,8))
    cb_plot=ax.scatter(df.fmelt,df.chamber,c=df.time,cmap='rainbow')
    ax_cbin =fig2.add_axes([0.13,0.78,0.23,0.03]) 
    cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
    ax_cbin.set_title('Myr')
    rrr1=fd.moving_window_smooth(df.ratio,10)
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
    print('-------plotting gravity-------')
    name='gravity_for_'+model
    fig, (ax,ax2)= plt.subplots(1,2,figsize=(22,12)) 
    dis,time,to,tom,fa,bg=get_gravity(1,end)
    qqq=ax.scatter(dis,time,c=fa,cmap='Spectral',vmax=400,vmin=-400)
    ax2.scatter(dis,time,c=bg,cmap='Spectral',vmax=400,vmin=-400)
    fig.colorbar(qqq,ax=ax)
    ax2.set_title('bourger gravoty anomaly')
    ax.set_title('free-air gravity anomaly')
    plt.savefig(figpath+model+'_gravity.png')
    print('=========== DONE =============')
if phase_plot:
    print('-----plotting single layer-----')
    name = 'phase_for'+model
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    colors = ["#CECCD0","#FF00FF","#8BFF8B","#7158FF","#FF966F",
          "#9F0042","#660000","#524B52","#D14309","#5AB245",
          "#004B00","#008B00","#455E45","#B89FCE","#C97BEA",
          "#525252","#FF0000","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    #rainbow = cm.get_cmap('Set1',20)
    #colors = rainbow(np.linspace(0, 1, 20))
    #phase15= matplotlib.colors.ListedColormap(colors)
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
    print('=========== DONE =============')
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
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_acc.png')
    print('=========== DONE =============')
if melting_plot:
    print('---plotting melting location---')
    name='melting_'+model
    #df=pd.read_csv(savepath+name+'.csv')
    time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    fig, (ax) = plt.subplots(1,1,figsize=(18,12))
    #ax.bar(time,phase_p9,width=0.17,color='orange',label='serpentinite ')
    ax.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='olivine')
    ax.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediments')
    ax.bar(time,phase_p3,bottom=phase_p10+phase_p4+phase_p10,width=0.17,color='#2360fa',label='basalt')
    ax.set_xlim(0,time[-1])
    ax.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('molten rocks (km3/km)',fontsize=20)
    ax.legend(fontsize=25)
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_bar_plot_melting.png')
    print('=========== DONE =============')
if force_plot_LR:
    print('-----plotting boundary force-----')
    filepath = savepath+model+'_forc.txt'
    fig, (ax,ax2)= plt.subplots(2,1,figsize=(12,8))   
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax.scatter(time,forc_l,c="#6A5ACD",s=4)
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
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig.savefig(figpath+model+'_forc.png')
    print('=========== DONE =============')
if force_plot_RF:
    print('----- plotting ringforce-----')
    filepath =savepath+model+'_forc.txt' 
    fig2, (ax3)= plt.subplots(1,1,figsize=(10,8))   
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    ax3.scatter(time,ringforce,c="#483D8B",s=4)
    ax3.set_xlim(0,time[-1])
    ax3.tick_params(axis='x', labelsize=16)
    ax3.tick_params(axis='y', labelsize=16)
    ax3.grid()
    fig2.savefig(figpath+model+'_ringforc.png')
    print('=========== DONE =============')
if vel_plot:
    print('-----plotting model velocity-----')
    filepath = savepath+model+'_forc.txt'
    temp1=np.loadtxt(filepath)
    nloop,time,forc_l,forc_r,ringforce,vl,vr,lstime,limit_force = temp1.T
    fig3, (ax4)= plt.subplots(1,1,figsize=(10,8))
    movvl = fd.moving_window_smooth(vl,500)
    #ax4.plot(time,vl*31545741325,c="darkred",lw=2)
    ax4.plot(time,movvl*31545741325,c="#000080",lw=2)
    ax4.set_xlim(0,time[-1])
    ax4.set_title('oceanic side velocity',fontsize=16)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.grid()
    ax4.set_xlabel('Time (Myr)',fontsize=16)
    ax4.set_ylabel('Velocity (mm/yr)',fontsize=16)
    bwith = 3
    ax4.spines['bottom'].set_linewidth(bwith)
    ax4.spines['top'].set_linewidth(bwith)
    ax4.spines['right'].set_linewidth(bwith)
    ax4.spines['left'].set_linewidth(bwith)
    fig3.savefig(figpath+model+'_vel.png')
    print('=========== DONE =============')
if stack_topo_plot:
    print('-----plotting topography-----')
    name=model+'_stack_topography.txt'
    xmean,ztop=np.loadtxt(savepath+name).T
    fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
    ax2.plot(xmean,ztop,c="#000080",lw=3)
    ax2.set_xlim(0,max(xmean)+10)
    bwith = 3
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig2.savefig(figpath+model+'_topo_analysis.png')
    print('=========== DONE =============')
if stack_gem_plot:
    print('-----plotting stacked geometry-----')
    name=model+'_stack_slab.txt'
    xmean,ztop=np.loadtxt(savepath+name).T
    fig2, (ax2) = plt.subplots(1,1,figsize=(8,6))
    ax2.plot(xmean,ztop,c="#000080",lw=3)
    ax2.set_xlim(0,max(xmean)+10)
    bwith = 3
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig2.savefig(figpath+model+'_gem.png')
    print('=========== DONE =============')
if flat_slab_plot:
    print('-----plotting flatslab-----')
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    fig2, (ax1,ax2) = plt.subplots(2,1,figsize=(10,8))
    ax1.plot(time,length,c="#000080",lw=3)
    ax2.plot(time,depth,c="#000080",lw=3)
    ax1.set_xlim(0,45)
    ax1.set_title('flat slab properties',fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.grid()
    ax1.set_ylabel('length (km)',fontsize=16)
    ax2.set_xlim(0,45)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.grid()
    ax2.set_xlabel('Time (Myr)',fontsize=16)
    ax2.set_ylabel('depth (km) ',fontsize=16)
    bwith = 3
    ax1.spines['bottom'].set_linewidth(bwith)
    ax1.spines['top'].set_linewidth(bwith)
    ax1.spines['right'].set_linewidth(bwith)
    ax1.spines['left'].set_linewidth(bwith)
    ax2.spines['bottom'].set_linewidth(bwith)
    ax2.spines['top'].set_linewidth(bwith)
    ax2.spines['right'].set_linewidth(bwith)
    ax2.spines['left'].set_linewidth(bwith)
    fig2.savefig(figpath+model+'_flatslab_length.png')
    print('=========== DONE =============')
