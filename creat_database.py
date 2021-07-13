#!/usr/bin/env python


import flac
import os,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import function_savedata as fs
import function_for_flac as f2
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


    ############## DO WHAT ###############
## creat data
trench = 0
magma = 0

# plot data
trench_plot = 0
magma_plot = 0
marker_number = 16 

#-------------------------------- setting ----------------------------------
path = '/home/jiching/geoflac/'
sys.path.append("/home/jiching/geoflac/util")
model = sys.argv[1]
os.chdir(path+model)

fl = flac.Flac()
end = fl.nrec
nex = fl.nx - 1
nez = fl.nz - 1
#----------------------------------------------------------------------------
def read_time(start,end_frame):
    timestep=[0]
    for step in range(start,end_frame+1):
        timestep.append(fl.time[step])
    # timestep=np.array(timestep)
    return timestep
time=read_time(1,end-1)
def trench(start=1,end_frame=end):
    trench_x=[0]
    trench_z=[0]
    for i in range(start,end_frame):
        x,z = fl.read_mesh(i+1)
        sx,sz=f2.get_topo(x,z,i+1)
        arc_ind,trench_ind=f2.find_trench_index(z)
        trench_x.append(sx[trench_ind])
        trench_z.append(sx[trench_ind])
    return trench_x,trench_z
trenchfile=path+'data/trench_for_'+model+'.csv'
def get_magma(start=1,end_frame=end-1):
    tol_melt=np.zeros(end_frame)
    tol_chamber=np.zeros(end_frame)
    for i in range(start,end_frame):
        x,z=fl.read_mesh(i)
        mm=fl.read_fmelt(i)
        cc=fl.read_chamber(i)
        mx,mz,mnum=f2.melt_element(x,z,i,mm)
        cx,cz,cnum=f2.chamber_element(x,z,i,cc)
        tol_melt[i]=0
        tol_chamber[i]=0
        for inn in range(len(mx)):
            x,z=fl.read_mesh(end)
            marea=f2.read_area(x,z,mx[inn],mz[inn])
            carea=f2.read_area(x,z,cx[inn],cz[inn])
            meltvol=marea * mnum[inn]
            chambervol=marea * cnum[inn]
            tol_melt[i] +=meltvol
            tol_chamber[i] +=chambervol
    return tol_melt,tol_chamber
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

#----------------------------------------------------------------------------
# if not os.path.exists(trenchfile):
if trench:
    print('-----creat trench database-----')
    name='trench_for_'+model
    trench_x,trench_z=trench()
    fs.save_3array(name,path+model,time,trench_x,trench_z,
                'time','trench_x','trench_z')
    print('=========== DONE =============')
# if not os.path.exists(magmafile):
if magma:
    print('-----creat magma database-----')
    name='magma_for_'+model
    tol_melt,tol_chamber=get_magma()      
    fs.save_3array(name,path+model,time,tol_melt,tol_chamber,
                   'time','fmelt','chamber')
    print('=========== DONE =============')

##-------------------------- plot --------------------------
if trench_plot:
    name='trench_for_'+model
    df = pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    ax.plot(df.trench_x,df.time,c='k',lw=2)
    plt.savefig(path+'figure/'+model+'_trench.jpg')
if magma_plot:
    name='magma_for_'+model
    df = pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax,ax2) = plt.subplots(2,1,figsize=(10,12))
    ax.bar(df.time,df.fmelt,width=0.3,color='tomato',label='fmelt')
    ax2.bar(df.time,df.chamber,width=0.3,color='orange',label='magma')
    plt.savefig(path+'figure/'+model+'_magma.jpg')
if marker_number != 0:
    mr = count_marker(marker_number)
    plt.plot(mr,c='b')
    