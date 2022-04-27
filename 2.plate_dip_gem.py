#!/usr/bin/env python
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as fd
from Main_creat_database import oceanic_slab,nodes_to_elements
model = str(sys.argv[1])
#model = 'Ref'
path = '/home/jiching/geoflac/'+model+'/'
#path = '/scratch2/jiching/03model/'+model+'/'
#path = '/scratch2/jiching/'+model+'/'
#path = 'F:/model/'+model+'/'
#path = 'D:/model/'+model+'/'
savepath = '/home/jiching/geoflac/data'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
fig, (ax)= plt.subplots(1,1,figsize=(17,12))

plot_dip = 0

phase_oceanic = 3
phase_ecolgite = 13
phase_oceanic_1 = 17
phase_ecolgite_1 = 18
angle = np.zeros(end)

rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))
stslab = 0;xmean=0;ictime=20;width=700
for i in range(1,end):
    crust_x,crust_z = oceanic_slab(i)
    x, z = fl.read_mesh(i)
    ele_x, ele_z = nodes_to_elements(x,z)
    phase = fl.read_phase(i)
    trench_ind = np.argmin(z[:,0]) 
    if z[:,0][trench_ind]>-2:
        print(i)
        continue
    x_trench = ele_x[:,0][np.argmin(ele_z[:,0])]
    within_plot = (ele_x[:,0]>x_trench)* (crust_z < 0)
    ax.plot(crust_x[within_plot]-(x_trench),crust_z[within_plot],color=newcolors[i],zorder=1)
    if not True in (crust_z < -80):
        continue
    if i >=end-ictime:
        stslab += crust_z
        xmean += (crust_x-x_trench)
    if plot_dip:
        ind_within_80km = (crust_z >= -80) * (crust_z < -5)
        crust_xmin = np.amin(crust_x[ind_within_80km])
        crust_xmax = np.amax(crust_x[ind_within_80km])
        crust_zmin = np.amin(crust_z[ind_within_80km])
        crust_zmax = np.amax(crust_z[ind_within_80km])
        dx = crust_xmax - crust_xmin
        dz = crust_zmax - crust_zmin
        angle[i] = math.degrees(math.atan(dz/dx))
xx=(xmean[within_plot]/ictime)
zz=(stslab[within_plot]/ictime)
ax.plot(xx[xx>0][:-1],zz[xx>0][:-1],c='green',lw=4)
ax.set_aspect('equal')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Depth (km)')
ax.set_xlim(-10,width)
ax.set_title('Geometry of Subducted Slab')
fs.save_2txt(str(model)+'_stack_slab',savepath,xx[xx>0],zz[xx>0])
fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_gem.jpg')

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
def oceanic_slab2(frame):
    bet = 2
    mx, mz, age, phase, ID, a1, a2, ntriag= fl.read_markers(frame)  
    trench_ind = np.argmin(z[:,0]) 
    x_trench,z_trench = x[trench_ind,0], z[trench_ind,0]
    x_ocean = mx[(phase==phase_ecolgite)+(phase==phase_oceanic)]
    z_ocean = mz[(phase==phase_ecolgite)+(phase==phase_oceanic)]
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
    return ox,oz

fff,aaa = plt.subplots(1,1,figsize=(17,12))
crust_x,crust_z = oceanic_slab(i)
within_plot = (ele_x[:,0]>x_trench)* (crust_z < 0)*(crust_z > -160)
xx = (crust_x-x_trench)[within_plot]
xx = fd.moving_window_smooth(xx,10)
zz = crust_z[within_plot]
zz = fd.moving_window_smooth(zz,10)
aaa.plot(xx,zz,color='#000080',lw=5)
aaa.set_aspect('equal')
aaa.set_xlabel('Distance (km)')
aaa.set_ylabel('Depth (km)')
aaa.set_xlim(0,width)
aaa.set_ylim(-150,0)
aaa.set_title('Geometry of Subducted Slab')
cx,cz = oceanic_slab2(i)
within_plot = (cz < 0)*(cz > -160)
xx = (cx-x_trench)
xx = fd.moving_window_smooth(xx,10)
zz = cz
zz = fd.moving_window_smooth(zz,10)
with_save = (xx<100)*(xx>0)
aaa.plot(xx[with_save],zz[with_save],color='tomato',lw=5)
aaa.set_aspect('equal')
fff.savefig('/home/jiching/geoflac/'+'figure/'+model+'_finalgem.jpg')
#fs.save_2txt(str(model)+'_final_slab',savepath,xx[with_save],zz[with_save])
if plot_dip:
    fig2, (ax2)= plt.subplots(1,1,figsize=(17,12))
    ax2.plot(fl.time[angle>0],angle[angle>0],c='blue',lw=2)
    ax2.set_xlim(0,fl.time[-1])
    ax2.set_title('Angle Variation')
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Angel ($^\circ$)')
    fig2.savefig('/home/jiching/geoflac/'+'figure/'+model+'_dip.jpg')
