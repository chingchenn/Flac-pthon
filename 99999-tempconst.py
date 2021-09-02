#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs



model = 'w0823'
#path = '/home/jiching/geoflac/'+model+'/'
path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
savepath='/Users/ji-chingchen/Desktop/model/data/'
print(model)
os.chdir(path)
# fl = flac.Flac()
fl = flac.FlacFromVTK();end = fl.nrec
melting_plot = 1
melting = 1


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
time =fl.time
def melting_phase(start_vts=1,model_steps=end-1):
    melt_num = np.zeros(end)
    phase_p4=np.zeros(end)
    phase_p9=np.zeros(end)
    phase_p10=np.zeros(end)
    po=np.zeros(end)
    for i in range(1,end):
        c=0;p9=0;p4=0;pk=0;p10=0
        x, z = fl.read_mesh(i)
        mm=fl.read_fmelt(i)
        phase=fl.read_phase(i)
        for xx in range(len(mm)):
            for zz in range(len(mm[0])):
                if mm[xx,zz] != 0:
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
if melting:
    print('-----creat magma database-----')
    name='melting_'+model
    phase_p4,phase_p9,phase_p10,po=melting_phase()
    fs.save_5array(name,savepath,time,phase_p4,phase_p9,phase_p10,po,
                'time','phase_4','phase_9','phase_10','others')
    print('=========== DONE =============')       
if melting_plot:
    path = '/Users/ji-chingchen/Desktop/model/'
    name='melting_'+model
    df=pd.read_csv(path+'data/'+name+'.csv')
    fig, (ax) = plt.subplots(1,1,figsize=(18,12))
    ax.bar(df.time,df.phase_9,width=0.17,color='orange',label='serpentinite ')
    ax.bar(df.time,df.phase_4,bottom=df.phase_9,width=0.17,color='seagreen',label='olivine')
    ax.bar(df.time,df.phase_10,bottom=df.phase_9+df.phase_4,width=0.17,color='tomato',label='sediments')
    ax.bar(df.time,df.others,bottom=df.phase_9+df.phase_4+df.phase_10,width=0.17,color='k',label='serpentinite')
    ax.set_xlim(0,24)
    ax.grid()
    ax.tick_params(axis='x', labelsize=16 )
    ax.tick_params(axis='y', labelsize=16 )
    ax.set_title('Model : '+model,fontsize=25)
    ax.set_xlabel('Time (Myr)',fontsize=20)
    ax.set_ylabel('number of elements in phases',fontsize=20)