#!/usr/bin/env python
import os
import flac
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

model_list=['w0940','w1013']
# name_list=['1.2e-12','6e-13 ','9e-13']
name_list=['1.6e-11','3.2e-11']
newcolors = ['#2F4F4F','#4682B4','#CD5C5C','#708090','#AE6378','#282130','#7E9680','#24788F','#849DAB','#EA5E51','#35838D','#4198B9','#414F67','#97795D','#6B0D47','#A80359','#52254F']
fig, (ax,ax2,ax3,ax4,ax5) = plt.subplots(5,1,figsize=(15,17))
for qq,model in enumerate(model_list):
    path = '/scratch2/jiching/'+model+'/'
    path = '/home/jiching/geoflac/'+model+'/'
    path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
    os.chdir(path)
    # fl = flac.Flac();end = fl.nrec
    fl = flac.FlacFromVTK();end = fl.nrec
    nex = fl.nx - 1;nez = fl.nz - 1
    print(model)    
    
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
    melt,magma,yymelt,yychamber,arc_vol=get_magma(1,end)
    
    ax.plot(time,yymelt,color=newcolors[qq],label=name_list[qq])
    ax2.plot(time,yychamber,color=newcolors[qq],label=name_list[qq])
    ax3.plot(time,melt,color=newcolors[qq],label=name_list[qq])
    ax4.plot(time,magma,color=newcolors[qq],label=name_list[qq])
    ax5.plot(time,arc_vol,color=newcolors[qq],label=name_list[qq])

ax5.set_xlabel('Time (Myr)',fontsize=20)
ax.set_ylabel('melt * area',fontsize=20)
ax2.set_ylabel('chamber *area',fontsize=20)
ax3.set_ylabel('max fmelt',fontsize=20)
ax4.set_ylabel('max fmagma',fontsize=20)
ax5.set_ylabel('arc area',fontsize=20)
ax.set_ylim(0,5)
ax2.set_ylim(0,100)
ax3.set_ylim(0,0.04)
ax4.set_ylim(0,0.051)
ax5.set_ylim(0,300)
ax.set_xlim(0,24)
ax2.set_xlim(0,24)
ax3.set_xlim(0,24)
ax4.set_xlim(0,24)
ax5.set_xlim(0,24)
ax.tick_params(axis='x', labelsize=16 )
ax2.tick_params(axis='x', labelsize=16 )
ax3.tick_params(axis='x', labelsize=16 )
ax4.tick_params(axis='x', labelsize=16 )
ax5.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16 )
ax2.tick_params(axis='y', labelsize=16 )
ax3.tick_params(axis='y', labelsize=16 )
ax4.tick_params(axis='y', labelsize=16 )
ax5.tick_params(axis='y', labelsize=16 )
ax5.legend(title = "productivity",fontsize = 16,loc=2)
ax.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
