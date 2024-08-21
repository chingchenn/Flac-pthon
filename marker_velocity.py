#!/usr/bin/env python
import flac
import sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_for_flac as f2
import function_for_flac as fd

# model = str(sys.argv[1])
model = 'Cocos_x0106'
# path = '/home/jiching/geoflac/'+model
#sys.path.append('/home/jiching/geoflac/util/')
path = '/Users/chingchen/Desktop/model/'
savepath = '/Users/chingchen/Desktop/data/'
figpath = '/Users/chingchen/Desktop/FLAC_Works/Observation_Mexico/'
os.chdir(path+model)
fl = flac.Flac();end = fl.nrec
plt.rcParams['figure.figsize'] =10,8
# number_of_marker = 50
# xmax,xmin=760,450
# zmax,zmin=-20,-50
# dis = np.zeros(end)
# t=np.linspace(1,end,end)
labelsize=20
fontsize=20
bwith=2
def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
# fig, (ax)= plt.subplots(1,1,figsize=(10,8))



fig3, (ax1)= plt.subplots(1,1,figsize=(15,7))  
#--------------------------------- FIG1 melting phase -------------------------
name='melting_'+model
# time,phase_p3,phase_p4,phase_p13,phase_p10 = np.loadtxt(savepath+name+'.txt').T
phase_p3=np.zeros(end)  # basalt
phase_p4=np.zeros(end)  # perditote
phase_p13=np.zeros(end)  # eclogite 
phase_p10=np.zeros(end)
for frame in range(1,end):
    c=0;p13=0;p4=0;p10=0;p3=0
    area = fl.read_area(frame)
    sed,bs,prid = fl.read_phase_melt(frame)
    mm=fl.read_fmelt(frame)
    x, z = fl.read_mesh(frame)
    ele_x, ele_z = nodes_to_elements(x,z) 
    # print(frame,np.max(sed),np.max(bs),np.max(prid))
    for xx in range(len(bs)):
        for zz in range(len(bs[0])-1):
            if mm[xx,zz] != 0:
                
                if sed[xx,zz] !=0:
                    p10 += area[xx,zz]*sed[xx,zz]/1e6
                if bs[xx,zz] !=0:
                    p3  += area[xx,zz]*bs[xx,zz]/1e6
                if prid[xx,zz] !=0:
                    p4  += area[xx,zz]*prid[xx,zz]/1e6
            # if mm[xx,zz] == 0 and prid[xx,zz] !=0:
            #     print(frame,ele_x[xx,zz],ele_z[xx,zz])
                
    phase_p3[frame]=p3
    phase_p4[frame]=p4
    phase_p10[frame]=p10
                
                    
# phase_p10 = sed
# phase_p3 = bs
# phase_p4 = prid
time = fl.time
total = (phase_p10+phase_p3+phase_p4)
ax1.bar(time,phase_p4,width=0.17,color='seagreen',label='peridotite')
ax1.bar(time,phase_p10,bottom=phase_p4,width=0.17,color='#B22222',label='sediment')
ax1.bar(time,phase_p3,bottom=phase_p4+phase_p10,width=0.17,color='#4169E1',label='eclogite')


axx=ax1.twinx()
kkk = fd.moving_window_smooth(phase_p3[total>0]/(phase_p10+phase_p3+phase_p4)[total>0]*100, 5)
kkk2 = fd.moving_window_smooth(phase_p4[total>0]/(phase_p10+phase_p3+phase_p4)[total>0]*100, 5)
axx.plot(time[total>0],kkk,c='#4169E1',lw=3)
axx.plot(time[total>0],kkk2,c='seagreen',lw=3)

for aaa in [ax1]:
    aaa.tick_params(labelsize=fontsize)
    aaa.grid()
    aaa.set_xlim(0,40)
    #aaa.axvspan(time[0],time[-1],facecolor='0.5', alpha=0.1)
    #aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
    for axis in ['top','bottom','left','right']:
        aaa.spines[axis].set_linewidth(bwith)
#------------------------------figure setting---------------------------------
# name='Nazca_a0702'+'_flatslab_time_len.txt'
# name=model+'_flatslab_time_len.txt'
# time_flat,length,depth=np.loadtxt(savepath+name).T
# axx.tick_params(labelsize=fontsize-5,labelcolor="#4169E1")
# axx.set_ylim(-20,100)
# axx.set_ylabel('% of melting rocks',fontsize=fontsize,color ="#4169E1")
# for aaa in [ax1]:
#     aaa.tick_params(labelsize=fontsize)
#     aaa.grid()
#     aaa.set_xlim(5,40)
#     # aaa.vlines(x=time_flat[0], ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
#     aaa.vlines(x=16.4, ymin=0, ymax=600, colors='purple', ls='--', lw=4,)
#     for axis in ['top','bottom','left','right']:
#         aaa.spines[axis].set_linewidth(bwith)
# ax1.set_ylim(0,40)    
# ax1.set_ylabel('molten rocks (km$^3$/km)',fontsize=fontsize)
# ax1.set_xlabel('time (Myr)',fontsize=fontsize)
# fig3.tight_layout()
# ax1.legend(fontsize=fontsize,facecolor='white')




# melt_num = np.zeros(end)
# phase_p3=np.zeros(end)  # basalt
# phase_p4=np.zeros(end)  # perditote
# phase_p13=np.zeros(end)  # eclogite 
# phase_p10=np.zeros(end) # sediment
# list_melt_marker_eclogite=[]
# for i in range(80,end):
#     mm=fl.read_fmelt(i)
#     phase=fl.read_phase(i)
#     area = fl.read_area(i)
#     x,z = fl.read_mesh(i)
#     ele_x, ele_z = nodes_to_elements(x,z)
#     for xx in range(len(mm)):
#         for zz in range(len(mm[0])-1):
#             if mm[xx,zz] != 0: 
            
#                 xp,zp,agep,phase_p ,ID_p, a1, a2, ntriag = fl.read_markers(i)
#                 x1 = x[xx,zz]
#                 x2 = x[xx+1,zz]
#                 x3 = x[xx,zz+1]
#                 x4 = x[xx+1,zz+1]
#                 z1 = z[xx,zz]
#                 z2 = z[xx+1,zz]
#                 z3 = z[xx,zz+1]
#                 z4 = z[xx+1,zz+1]
                
#                 xmin=min(x1,x2,x3,x4)
#                 xmax=max(x1,x2,x3,x4)
#                 zmin=min(z1,z2,z3,z4)
#                 zmax=max(z1,z2,z3,z4)
#                 melt_marker_eclogite = ID_p[(xp>=xmin)*(xp>=xmax)*(zp<z3)*(zp>=z4)*(phase_p==13)]
#                 if len(melt_marker_eclogite)!=0:
#                     for yy in range(len(melt_marker_eclogite)):
#                         qq = melt_marker_eclogite[yy]
#                         # print(qq)
#                         list_melt_marker_eclogite.append(qq)
                
# fig2,ax2 = plt.subplots(1,1,figsize=(10,10))
# x = np.linspace(0,514)
# y = -0.0375 * x + 20.1
# basalt_change = (50/255, 200/255, 180/255)
# ax2.plot(x,y,c=basalt_change,lw=5)
# x = np.linspace(515,1300)
# y = 0.0022 * x - 0.3
# ax2.plot(x,y,c=basalt_change,lw=5,label='basalt-eclogite')
# pressure_limit = 5 # GPa
# pressure=np.linspace(0,pressure_limit,100)
# sss=np.zeros(len(pressure))
# for q,dd in enumerate(pressure):
#     if dd<1:
#         ss=1050-420*(1-np.exp(-dd*3.3))
#     elif dd>2.38:
#         ss=(dd+14)*43
#     else:
#         ss=630+26*((-dd)**2)/2 
#     sss[q] = ss
# ax2.plot(sss,pressure,c='#FF9900',lw=5,label='solidus')

# x = np.linspace(710,1050)
# y = -1.25/350*x+5
# ax2.plot(x,y,c='#FF9900',lw=5) # solidus

# x = np.linspace(680,1050)
# y = 0.65/400*x-0.45625
# ax2.plot(x,y,c='#FF9900',lw=5) # solidus

# ax2.set_ylim(0,pressure_limit)
# ax2.set_xlim(0,1200)
# axdep = ax2.twinx()
# axdep.set_ylim(0,pressure_limit*1e9/3300/10/1e3)
# ax2.legend(fontsize=labelsize)
# axdep.tick_params(axis='y', labelsize=labelsize)
# ax2.tick_params(labelsize=labelsize)
# ax2.set_xlabel('Temperature ($^\circ$C)',fontsize=labelsize)
# ax2.set_ylabel('Pressure (GPa)',fontsize=labelsize)
# axdep.set_ylabel('Depth (km)',fontsize=labelsize)
# for axis in ['top','bottom','left','right']:
#     ax2.spines[axis].set_linewidth(bwith)

            
# for kk,idd in enumerate(list_melt_marker_eclogite):
#     # px = np.zeros(end)
#     pz = np.zeros(end)
#     pPre = np.zeros(end)
#     pTemp = np.zeros(end)
#     pmark= np.zeros(end)
#     for ii,qq in enumerate(range(21,150)):
#         mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(qq)
#         temp = fl.read_temperature(qq)
#         if len(mx[ID == idd])==0:
#             continue
#         # px[ii] = mx[ID == idd]
#         pz[ii] = mz[ID == idd]
#         pmark[ii] = phase[ID == idd]
#         pPre[ii] = 3300*10*(-pz[ii])*1000/1e9
#         pTemp[ii] = flac.marker_interpolate_node(ntriag[ID==idd], a1[ID==idd], a2[ID==idd], fl.nz, temp)
#         pTemp[ii] = pTemp[ii] - 0.4*pz[ii]
#     # ax2.scatter(pTemp[pTemp>0],pPre[pTemp>0],c=pmark[pTemp>0],cmap = phase15,s = 40,vmin = 1,vmax = 20)
#     # axdep.scatter(pTemp[pTemp>0],-pz[pTemp>0],c=pmark[pTemp>0],cmap = phase15,s = 40,vmin = 1,vmax = 20)
#     ax2.plot(pTemp[pTemp>0]+50,pPre[pTemp>0],c='purple')                  
    
