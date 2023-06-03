#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as fd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#=========================setting=============================
# model = str(sys.argv[1])
# model = 'Nazca_0502'
# path = '/home/jiching/geoflac/'+model+'/'
# path = '/Users/chingchen/Desktop/model/'+model+'/'
path = '/Users/chingchen/Desktop/model/'

fig5 = 0
fig6 = 0
fig7 = 1
fig8 = 0

plt.rcParams["font.family"] = "Times New Roman"
# model='Ref_Nazca'
model='Ref_Cocos'
# label_list=['110 km','120 km','130 km','140 km','150 km','160 km']
newcolors = ['#2F4F4F','#A80359','#4198B9','#AE6378',
             '#35838D','#97795D','#7E9680','#4682B4',
             '#708090','#282130','#24788F','#849DAB',
             '#EA5E51','#414F67','#6B0D47','#52254F'] 
savepath='/home/jiching/geoflac/data/'
savepath = '/Users/chingchen/Desktop/data/'

#========================= Function =========================
def get_magma(end):
#=========================Time Series=========================
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    rrr=np.zeros(end)
    depth_magma = np.zeros(end)
#=========================main code===========================
    for i in range(1,end):
        x,z = fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        chamber=fl.read_fmagma(i)
        melt[i] = np.max(mm)
        magma[i] = np.max(chamber)
        ele_x,ele_z = flac.elem_coord(x, z)
        ui,uj= np.unravel_index(chamber.argmax(), chamber.shape)
        depth_magma[i] = -ele_z[ui,uj]
        if  magma[i]!=0:
            rrr[i]= melt[i]/magma[i]
        arc_vol[i]=np.sum(fl.read_area(i)[phase ==14])/1e6
        yymelt[i]=(fl.read_fmelt(i)*fl.read_area(i)/1e6).sum()
        yychamber[i]=(fl.read_fmagma(i)*fl.read_area(i)/1e6).sum()
    return melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma

#=============================================================
###=======================plot================================
#=============================================================
# fig, (ax,ax2,ax3,ax4,ax5)= plt.subplots(5,1,figsize=(25,18))
# ax.plot(fl.time,yymelt,color='tomato')
# ax2.plot(fl.time,yychamber,color='orange')
# ax3.bar(fl.time,melt,width=0.1,color='tomato')
# ax4.bar(fl.time,magma,width=0.1,color='orange')
# ax5.plot(fl.time,arc_vol,color='orange',label='magma')
# ax.set_title('Model : '+model,fontsize=25)
# ax.set_ylabel('melt * area',fontsize=20)
# ax2.set_ylabel('magma friction * area',fontsize=20)
# ax3.set_ylabel('max melt fraction',fontsize=20)
# ax4.set_ylabel('max chamber fraction',fontsize=20)
# ax5.set_xlabel('Time (Myr)',fontsize=20)
# ax5.set_ylabel('arc area',fontsize=20)
# #ax.set_ylim(0,0.8)
# #ax2.set_ylim(0,10*1e-3)
# #ax3.set_ylim(0,10*1e-3)
# #ax4.set_ylim(0,3*1e-5)
# #ax5.set_ylim(0,300)
# ax.set_xlim(0,fl.time[-1])
# ax2.set_xlim(0,fl.time[-1])
# ax3.set_xlim(0,fl.time[-1])
# ax4.set_xlim(0,fl.time[-1])
# ax.grid()
# ax2.grid()
# ax3.grid()
# ax4.grid()
# ax5.grid()
# ax.tick_params(axis='x', labelsize=16)
# ax2.tick_params(axis='x', labelsize=16)
# ax3.tick_params(axis='x', labelsize=16)
# ax4.tick_params(axis='x', labelsize=16)
# ax5.tick_params(axis='x', labelsize=16 )
# ax.tick_params(axis='y', labelsize=16)
# ax2.tick_params(axis='y', labelsize=16)
# ax3.tick_params(axis='y', labelsize=16)
# ax4.tick_params(axis='y', labelsize=16)
# ax5.tick_params(axis='y', labelsize=16)
# #-------------------------------------------------------------
# fig2,(ax,ax2)=plt.subplots(1,2,figsize=(25,18))
# cb_plot=ax.scatter(melt,magma,c=fl.time,cmap='rainbow')
# ax_cbin = fig2.add_axes([0.13, 0.78, 0.23, 0.03])
# cb = fig2.colorbar(cb_plot,cax=ax_cbin,orientation='horizontal')
# ax_cbin.set_title('Myr (km)')
# rrr1=f2.moving_window_smooth(rrr,12)
# ax2.plot(fl.time,rrr1,color='k',lw=3)
# ax2.plot(fl.time,rrr,color='gray',linestyle=':')
# fig.savefig('/home/jiching/geoflac/figure/'+model+'_magma_parameter_time_series.png')
# fig2.savefig('/home/jiching/geoflac/figure/'+model+'_max_ratio.png')

bwith=3


fontsize = 30
if fig5:
    fig5, (ax2,ax1,ax3)= plt.subplots(3,1,figsize=(12,16),gridspec_kw={'height_ratios':[1,1,1]})  
    model='Nazca_a0702'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-2]
        hhh = (melting>1e-2)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=10,cmap='OrRd',vmax=3,vmin=-0)

    melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma=get_magma(end)
    time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    qqq=ax2.scatter(time[melt>1e-3],xmelt[melt>1e-3],c=np.log10(melt[melt>1e-3]*100),s=50,cmap='OrRd',vmax=2,vmin=-2)
    qqq=ax2.scatter(time[melt>1e-3],xmelt[melt>1e-3],c=melt[melt>1e-3]*100,s=50,cmap='OrRd',vmax=3,vmin=0)
    divider = make_axes_locatable(ax2)
    cax = plt.axes([0.99, 0.7, 0.02, 0.28])
    cbar=fig5.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('Melting Percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p9,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    ax1.bar(time,phase_p4+phase_p9,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediment')
    ax1.bar(time,phase_p3*0.1,bottom=phase_p4+phase_p9+phase_p10,width=0.17,color='darkblue',label='eclogite')
    
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    qqq1=ax3.scatter(time[yychamber>1e-3],yychamber[yychamber>1e-3],c=depth_magma[yychamber>1e-3],cmap='gist_earth_r',vmin=70,vmax=90)
    divider = make_axes_locatable(ax2)
    cax = plt.axes([0.99, 0.06, 0.02, 0.28])
    cbar=fig5.colorbar(qqq1,ax=ax3,cax=cax,orientation='vertical')
    cbar.set_label('Chamber Depth (km)',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    
    #================================figure setting================================
    # ax2.set_title('Partial melting distance from trench v.s. time ' ,fontsize=fontsize)
    # ax1.set_title('Molten rocks v.s. time ' ,fontsize=fontsize)
    # ax3.set_title('Chamber volume v.s. time ' ,fontsize=fontsize)
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1,ax3]:
        aaa.tick_params(axis='x', labelsize=fontsize)
        aaa.tick_params(axis='y', labelsize=fontsize)
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.grid()
        aaa.set_xlim(0,40)
        aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        
    ax2.set_ylim(0,600)
    ax1.set_ylim(0,10)
    ax3.set_ylim(0,0.20)
    ax2.set_ylabel('Distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    ax3.set_ylabel('Chamber volume (km$^3$/km)',fontsize=fontsize)
    ax3.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig5.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    # fig5.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Ref_Nazca/melting_time_series_5.pdf')

if fig6:
    fig6, (ax2,ax1,ax3)= plt.subplots(3,1,figsize=(12,16),gridspec_kw={'height_ratios':[1,1,1]})  
    # fig7, (ax4,ax5)= plt.subplots(2,1,figsize=(12,8))  
    
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    time,trench_index, trench_x, trench_z = np.loadtxt(savepath+'trench_for_'+str(model)+'.txt').T
    for ii in range(1,end):
        melting = fl.read_fmelt(ii)
        x, z = fl.read_mesh(ii)
        ele_x, ele_z = flac.elem_coord(x,z)
        hod = ele_x[melting>1e-3]
        hhh = (melting>1e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=40,cmap='OrRd',vmax=1,vmin=-0)
        hod = ele_x[melting>6e-3]
        hhh = (melting>6e-3)
        if len(hod)>0:
            ttt = np.ones(len(hod))*ii*0.2
            qqq=ax2.scatter(ttt, ele_x[hhh]-trench_x[ii-1],c =melting[hhh]*100 ,s=40,cmap='OrRd',vmax=1,vmin=-0)
    melt,magma,yymelt,yychamber,arc_vol,rrr,depth_magma=get_magma(end)
    time,melt,xmelt,zmelt=np.loadtxt(savepath+'metloc_for_'+model+'.txt').T
    # qqq=ax2.scatter(time[melt>1e-3],xmelt[melt>1e-3],c=melt[melt>1e-3]*100,s=50,cmap='OrRd',vmax=3,vmin=-0)
    # divider = make_axes_locatable(ax2)
    cax = plt.axes([0.99, 0.7, 0.02, 0.28])
    cbar=fig6.colorbar(qqq,ax=ax2,cax=cax,orientation='vertical')
    cbar.set_label('Melting percentages',fontsize=fontsize)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='tomato',label='sediment')
    # ax1.bar(time,phase_p3,bottom=phase_p4+phase_p10,width=0.17,color='#000080',label='basalt')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#000080',label='eclogite')
    
    ppptime = time
    
    
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    qqq1=ax3.scatter(time[yychamber>1e-3],yychamber[yychamber>1e-3],c=depth_magma[yychamber>1e-3],cmap='gist_earth_r',vmin=45,vmax=65)
    # qqq1=ax3.scatter(time[total>0],yychamber[total>0],c=kkk,cmap='brg',vmin=0,vmax=60)
    
    divider = make_axes_locatable(ax2)
    cax = plt.axes([0.99, 0.06, 0.02, 0.28])
    cbar=fig6.colorbar(qqq1,ax=ax3,cax=cax,orientation='vertical')
    cbar.set_label('Eclogite melting percentages',fontsize=fontsize+1)
    cbar.ax.tick_params(axis='y', labelsize=fontsize-2)
    cbar.ax.yaxis.set_label_position('right')
    
    #================================figure setting================================
    # ax2.set_title('Partial melting distance from trench v.s. time ' ,fontsize=fontsize)
    # ax1.set_title('Molten rocks v.s. time ' ,fontsize=fontsize)
    # ax3.set_title('Chamber volume v.s. time ' ,fontsize=fontsize)
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax2,ax1,ax3]:
        aaa.tick_params(axis='x', labelsize=fontsize)
        aaa.tick_params(axis='y', labelsize=fontsize)
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.grid()
        aaa.set_xlim(0,40)
        aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        
    ax2.set_ylim(0,300)
    ax3.set_ylim(0,20)
    # ax3.set_ylim(0,0.2)
    ax2.set_ylabel('Distance from trench (km)',fontsize=fontsize)
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    ax3.set_ylabel('Chambers volume (km$^3$/km)',fontsize=fontsize)
    ax3.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig6.tight_layout()
    
    
    # ax4.bar(ppptime,phase_p4+phase_p9,width=0.17,color='seagreen',label='peridotite')
    # ax4.bar(ppptime,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediment')
    # ax4.bar(ppptime,phase_p3,bottom=phase_p4+phase_p9+phase_p10,width=0.17,color='darkblue',label='basalt')
    # ax4.set_ylim(0,0.5)
    # kkk = fd.moving_window_smooth(phase_p3/(phase_p10+phase_p3+phase_p4)*100, 5)
    # ax5.bar(ppptime,kkk,color='darkblue',width=0.2)
    # # plt.close(fig7)
    ax1.legend(fontsize=fontsize,facecolor='white')
    # fig6.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Ref_Cocos/melting_time_series_v2.pdf')

if fig7:
    fig7, (ax1,ax3)= plt.subplots(2,1,figsize=(15,10))  
    # fig7, (ax4,ax5)= plt.subplots(2,1,figsize=(12,8))  
    
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='tomato',label='sediment')
    # ax1.bar(time,phase_p3,bottom=phase_p4+phase_p10,width=0.17,color='#000080',label='basalt')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#000080',label='eclogite')
    
    ppptime = time
    
    
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    axx=ax3.twinx()
    ax3.bar(time[total>0],yychamber[total>0],width=0.17,color='orange',alpha=0.5)
    axx.scatter(time[total>0],kkk,c='b')
    
    #================================figure setting================================
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax1,ax3,axx]:
        aaa.tick_params(axis='x', labelsize=fontsize)
        aaa.tick_params(axis='y', labelsize=fontsize)
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.set_xlim(13,40)
    ax1.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        
    
    ax1.set_ylim(0,0.75)
    ax3.set_ylim(0,7)
    axx.set_ylim(0,80)
    axx.grid()
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    ax3.set_ylabel('Chamber volume (km$^3$/km)',fontsize=fontsize)
    ax3.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig7.tight_layout()
    ax1.legend(fontsize=fontsize,facecolor='white')
    # fig7.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Discussion/slab_melting2.pdf')
    
if fig8:
    fig8, (ax1)= plt.subplots(1,1,figsize=(15,10))  
    # fig7, (ax4,ax5)= plt.subplots(2,1,figsize=(12,8))  
    
    model='Ref_Cocos'
    os.chdir(path+model)
    fl = flac.Flac();end = fl.nrec
    name='melting_'+model
    time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
    total = (phase_p10+phase_p3+phase_p4)
    kkk = fd.moving_window_smooth(phase_p13[total>0]/(phase_p10+phase_p13+phase_p4)[total>0]*100, 5)
    ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
    ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='tomato',label='sediment')
    # ax1.bar(time,phase_p3,bottom=phase_p4+phase_p10,width=0.17,color='#000080',label='basalt')
    ax1.bar(time,phase_p13,bottom=phase_p4+phase_p10+phase_p3,width=0.17,color='#000080',label='eclogite')
    
    ppptime = time
    
    
    name = 'magma_for_'+model+'.txt'
    temp1 = np.loadtxt(savepath+name)
    melt,chamber,yymelt,yychamber,rrr = temp1.T
    axx=ax1.twinx()
    # ax3.bar(time[total>0],yychamber[total>0],width=0.17,color='orange',alpha=0.5)
    axx.scatter(time[total>0],kkk,c='#8A2BE2')
    
    #================================figure setting================================
    # ax1.set_title('Molten rocks v.s. time ' ,fontsize=fontsize)
    # ax3.set_title('Chamber volume v.s. time ' ,fontsize=fontsize)
    name=model+'_flatslab_time_len.txt'
    time,length,depth=np.loadtxt(savepath+name).T
    for aaa in [ax1,axx]:
        aaa.tick_params(axis='x', labelsize=fontsize)
        aaa.tick_params(axis='y', labelsize=fontsize)
        aaa.spines['bottom'].set_linewidth(bwith)
        aaa.spines['top'].set_linewidth(bwith)
        aaa.spines['right'].set_linewidth(bwith)
        aaa.spines['left'].set_linewidth(bwith)
        aaa.set_xlim(13,40)
    ax1.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
        
    
    ax1.set_ylim(0,0.75)
    ax3.set_ylim(0,7)
    axx.set_ylim(0,80)
    axx.grid()
    ax1.set_ylabel('Molten rocks (km$^3$/km)',fontsize=fontsize)
    ax3.set_ylabel('Chamber volume (km$^3$/km)',fontsize=fontsize)
    ax3.set_xlabel('Time (Myr)',fontsize=fontsize)
    fig7.tight_layout()
    
    
    # ax4.bar(ppptime,phase_p4+phase_p9,width=0.17,color='seagreen',label='peridotite')
    # ax4.bar(ppptime,phase_p10,bottom=phase_p4+phase_p9,width=0.17,color='tomato',label='sediment')
    # ax4.bar(ppptime,phase_p3,bottom=phase_p4+phase_p9+phase_p10,width=0.17,color='darkblue',label='basalt')
    # ax4.set_ylim(0,0.5)
    # kkk = fd.moving_window_smooth(phase_p3/(phase_p10+phase_p3+phase_p4)*100, 5)
    # ax5.bar(ppptime,kkk,color='darkblue',width=0.2)
    # # plt.close(fig7)
    ax1.legend(fontsize=fontsize,facecolor='white')
    # fig7.savefig('/Users/chingchen/Library/CloudStorage/OneDrive-國立台灣大學/Thesis_figure/Discussion/slab_melting2.pdf')    
