#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as f2
model_list=['Chi01','chi0401','chi0402','chi0403','chi0404','chi0405']
name_list=['2','3','5','10','20','30']
rainbow = cm.get_cmap('rainbow',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
# depth1=-5
# depth2=-150
# case =1
if len(sys.argv) <= 1:
    print('''
          input [1] = plot style : 1 = plot dip only
                                  : 2 = plot dip with smoothing
          input [2] [3] = dip depth range from [2] to [3]
                          Note that [2] [3] is negative
          ''')
    sys.exit()
case =sys.argv[1]
depth1=sys.argv[2]
depth2=sys.argv[3]
if int(case)==1:
    print(11111111111111111)
    fig, (ax2)= plt.subplots(1,1,figsize=(13,8))
if int(case)==2: 
    print(22222222222222222)
    fig, (ax,ax2)= plt.subplots(2,1,figsize=(10,10))
for qq,model in enumerate(model_list):
    path = '/scratch2/jiching/22winter/'+model+'/'
    path = '/home/jiching/geoflac/'+model+'/'
    path = '/scratch2/jiching/03model/'+model+'/'
    #path = 'D:/model/'+model+'/'
    #path = '/scratch2/jiching/sem02model/'+model+'/'
    #path = '/scratch/jiching/summer2021/week11/'+model+'/'
    #path = '/scratch2/jiching/'+model+'/'
    # path = '/Volumes/My Book/model/'+model+'/'
    #path = '/Volumes/SSD500/model/'+model+'/'
    os.chdir(path)
    fl = flac.Flac();end = fl.nrec
    nex = fl.nx - 1;nez = fl.nz - 1
    
    phase_oceanic = 3
    phase_ecolgite = 13
    phase_oceanic_1 = 17
    phase_ecolgite_1 = 18
    angle = np.zeros(end)


    for i in range(1,end):
        x, z = fl.read_mesh(i)
        ele_x = (x[:fl.nx-1,:fl.nz-1] + x[1:,:fl.nz-1] + x[1:,1:] + x[:fl.nx-1,1:]) / 4.
        ele_z = (z[:fl.nx-1,:fl.nz-1] + z[1:,:fl.nz-1] + z[1:,1:] + z[:fl.nx-1,1:]) / 4.
        phase = fl.read_phase(i)
        trench_ind = np.argmin(z[:,0]) 
        crust_x = np.zeros(nex)
        crust_z = np.zeros(nex)
        for j in range(trench_ind,nex):
            ind_oceanic = (phase[j,:] == phase_oceanic) + (phase[j,:] == phase_ecolgite)+(phase[j,:] == phase_oceanic_1) + (phase[j,:] == phase_ecolgite_1)
            if True in ind_oceanic:
                crust_x[j] = np.average(ele_x[j,ind_oceanic])
                crust_z[j] = np.average(ele_z[j,ind_oceanic])

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
    if int(case)==2:
        nnewangle = f2.moving_window_smooth(angle[angle>0],3)
        ax.plot(fl.time[angle>0],nnewangle,c=newcolors[qq],lw=3,label=name_list[qq])
        # ax.legend(title = "convergent velocity mm/year")
        ax.set_xlim(0,fl.time[-1])
        ax.set_title('Smoothing Angle Variation')
        ax.set_ylabel('Angel ($^\circ$)')
    ax2.plot(fl.time[angle>0],angle[angle>0],c=newcolors[qq],lw=3,label=name_list[qq])
    ax2.legend(title = "geology zone",fontsize = 16)
    ax2.set_xlim(0,fl.time[-1])
    ax2.grid(axis='y')
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Angel ($^\circ$)')
    ax2.set_title('Angle Variation')
plt.savefig('/home/jiching/geoflac/figure/'+'dip_'+model_list[0]+'_'+model_list[1]+'_'+str(depth2)+'km'+'_.png')
