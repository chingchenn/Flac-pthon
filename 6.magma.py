#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as f2
import matplotlib.pyplot as plt
#=========================setting=============================
model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model+'/'
#path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
#=========================Time Series=========================
def get_magma(start_vts=1,model_steps=end-1):
    #=========================Time Series=========================
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    rrr=np.zeros(end)
    #=========================main code===========================
    for i in range(1,end):
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        chamber=fl.read_fmagma(i)
        melt[i] = np.max(mm)
        magma[i] = np.max(chamber)
        if  magma[i]!=0:
            rrr[i]= melt[i]/magma[i]
        arc_vol[i]=np.sum(fl.read_area(i)[phase ==14])/1e6
        yymelt[i]=(fl.read_fmelt(i)*fl.read_area(i)/1e6).sum()
        yychamber[i]=(fl.read_fmagma(i)*fl.read_area(i)/1e6).sum()
    return melt,magma,yymelt,yychamber,arc_vol,rrr
melt,magma,yymelt,yychamber,arc_vol,rrr=get_magma(1,end)
#=============================================================
###=======================plot================================
#=============================================================
fig, (ax,ax2,ax3,ax4,ax5)= plt.subplots(5,1,figsize=(25,18))
ax.plot(fl.time,yymelt,color='tomato')
ax2.plot(fl.time,yychamber,color='orange')
ax3.bar(fl.time,melt,width=0.1,color='tomato')
ax4.bar(fl.time,magma,width=0.1,color='orange')
ax5.plot(fl.time,arc_vol,color='orange',label='magma')
ax.set_title('Model : '+model,fontsize=25)
ax.set_ylabel('melt * area',fontsize=20)
ax2.set_ylabel('magma friction * area',fontsize=20)
ax3.set_ylabel('max melt fraction',fontsize=20)
ax4.set_ylabel('max chamber fraction',fontsize=20)
ax5.set_xlabel('Time (Myr)',fontsize=20)
ax5.set_ylabel('arc area',fontsize=20)
#ax.set_ylim(0,0.8)
#ax2.set_ylim(0,10*1e-3)
#ax3.set_ylim(0,10*1e-3)
#ax4.set_ylim(0,3*1e-5)
#ax5.set_ylim(0,300)
ax.set_xlim(0,fl.time[-1])
ax2.set_xlim(0,fl.time[-1])
ax3.set_xlim(0,fl.time[-1])
ax4.set_xlim(0,fl.time[-1])
ax.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='x', labelsize=16)
ax5.tick_params(axis='x', labelsize=16 )
ax.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax5.tick_params(axis='y', labelsize=16)
#-------------------------------------------------------------
fig2,(ax,ax2)=plt.subplots(1,2,figsize=(25,18))
cb_plot=ax.scatter(melt,magma,c=fl.time,cmap='rainbow')
ax_cbin = fig2.add_axes([0.13, 0.78, 0.23, 0.03])
cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
ax_cbin.set_title('Myr (km)')
rrr1=f2.moving_window_smooth(rrr,12)
ax2.plot(fl.time,rrr1,color='k',lw=3)
ax2.plot(fl.time,rrr,color='gray',linestyle=':')
fig.savefig('/home/jiching/geoflac/figure/'+model+'_magma_parameter_time_series.png')
fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')
