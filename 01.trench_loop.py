#!/usr/bin/env python
import flac
import os,sys
import numpy as np
import matplotlib.pyplot as plt
all_list=[['w0330','w0324']]

for yind,pair in enumerate(all_list): 
    model_list=all_list[yind]
    fig,ax1=plt.subplots(1,1,figsize=(12,8))
    newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
    for qq,model in enumerate(model_list):    
        path='/home/jiching/geoflac/'+model+'/'
        os.chdir(path)
        fl=flac.Flac()
        end=fl.nrec
        trench_x = np.zeros(end)
        trench_t = np.zeros(end)    
        for i in range(end):
            x,z = fl.read_mesh(i+1)
            phase = fl.read_phase(i+1)
            xtop=x[:,0]
            ztop=z[:,0]
            trench_id=np.argmin(ztop)
            if z[trench_id,0] > -2 :
                continue
            t=np.ones(xtop.shape)
            t[:]=i*0.2
            trench_t[i] = t[0]
            trench_x[i] = xtop[trench_id]
            # ax1.scatter(xtop,t,c=ztop,cmap=cmap,vmin=-20,vmax=10)
        ax1.plot(trench_x,trench_t,c=newcolors[qq],lw='4',label=str(model_list[qq]))
        ax1.legend()
        ax1.set_ylim(0,t[0])
        ax1.set_ylabel("Time (Ma)")
        ax1.set_xlabel("Distance (km)")
    plt.savefig('/home/jiching/geoflac/figure'+'/'+str(all_list[yind])+'_conpair_topo.jpg')
