#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:15:55 2022

@author: chingchen
"""

import sys, os
import numpy as np

import flac
import flac_interpolate as fi
# import flac_gravity3 as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

#-------------------------------------------------------------------
model = sys.argv[1]
frame = int(sys.argv[2])
# model = 'Nazca_a0624'
# frame = 10
plt.rcParams["font.family"] = "Times New Roman"
path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
# path = '/Users/chingchen/Desktop/model/'
savepath='/scratch2/jiching/data/'
#savepath='/Volumes/SSD500/data/'
# savepath='/home/jiching/geoflac/data/'
# savepath = '/Users/chingchen/Desktop/data/'
figpath='/scratch2/jiching/figure/'
# figpath = '/Users/chingchen/Desktop/figure/'

os.chdir(path+model)
fl = flac.Flac()
end = fl.nrec


phasein      = 1
cpp          = 1
grd          = 1
figure_plot2 = 0
import time
start_time = time.time()
# -------------------------------------------------------------------
# domain bounds
left = -300
right = 1000
up = 10
down = -300
dx = 1.2
dz = 1.1

def find_trench_index(z):
    '''Returns the i index of trench location.'''
    zz = z[:,0]
    # the highest point defines the forearc
    imax = zz.argmax()
    # the trench is the lowest point west of forearc
    i = zz[:imax].argmin()
    return i

def interpolate_phase(frame, xtrench):
    # domain bounds in km
    fi.xmin = xtrench + left
    fi.xmax = xtrench + right
    fi.zmin = down
    fi.zmax = up

    # resolution in km
    fi.dx = dx
    fi.dz = dz

    xx, zz, ph = fi.interpolate(frame, 'phase')
    return xx, zz, ph
frame_list=[26,51,76,101,126,150]
for frame in frame_list:
    x, z = fl.read_mesh(frame)
    itrench = find_trench_index(z)
    xtrench = x[itrench,0]
    ###############################################
    
    # get interpolated phase either from previous run or from original data
    phasefile = 'intp3-phase.%d' % frame
    phgrd = model+'_phase3.%d.grd' % frame
    if phasein:
        xx, zz, ph = interpolate_phase(frame, xtrench)
        f = open(phasefile, 'w')
        f.write('%d %d\n' % xx.shape)
        flac.printing(xx, zz, ph, stream=f)
        f.close()
    
    if grd:
        xmin = 250
        xmax = 1000
        zmin = down
        zmax = up
        
        cmd = 'tail -n +2 %(phasefile)s | gmt xyz2grd -G%(phgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f' % locals()
        print(cmd)
        os.system(cmd)
        
    ###############
    if cpp:
        cpcmd = ''' awk '{print $1,$2,$3+0}' %(phasefile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(phasefile)s.txt
        mv  %(phgrd)s %(savepath)s
        mv  %(model)s_%(phasefile)s.txt %(savepath)s
    ''' % locals()
        os.system(cpcmd)
    colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#DC143C","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#CD5C5C","#00FF00","#FFFF00","#7158FF"]
    phase15= matplotlib.colors.ListedColormap(colors)
    print('---------------end of '+str(frame)+'---------------')
    
    if figure_plot2:
        
        fig, (ax)= plt.subplots(1,1,figsize=(17,16))
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)
        bwith = 3
        #--------------------- phase plotting -------------------------
        from netCDF4 import Dataset
        # data = Dataset(savepath+model+'_phase3.'+str(frame)+'.grd', mode='r')
        # data = Dataset(savepath+model+'.grd', mode='r')
        data = Dataset(savepath+'NNN.grd', mode='r')
        x = data.variables['x'][:]
        z = data.variables['y'][:]
        ph = data.variables['z'][:]
        phh=ph.data[ph.data>0]
        xv, zv = np.meshgrid(x, z)
        ax.pcolormesh(xv,-zv,ph,cmap=phase15,vmin=1, vmax=20)
        # filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
        # x,z,ph=np.loadtxt(filepath).T
        # ax.scatter(x,-z,c = ph,cmap = phase15,vmax=19,vmin=1,s=1.5)
        # ax.scatter()
        ax.contour(xt,-zt,temp,cmap='rainbow',levels =[200,400,600,800,1000,1200],linewidths=3)
        # ---------------------- plot setting --------------------------
        ax.set_aspect('equal')
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.tick_params(axis='x', labelsize=23)
        ax.tick_params(axis='y', labelsize=23)
        ymajor_ticks = np.linspace(200,0,num=5)
        ax.set_yticks(ymajor_ticks)
        #xmajor_ticks = np.linspace(250,1000,num=6)
        #ax.set_xticks(xmajor_ticks)
        #ax.set_xlim(250,1000)
        ax.set_ylim(200,-30)
        xmajor_ticks = np.linspace(250,1000,num=7)
        ax.set_xticks(xmajor_ticks)
        ax.set_xlim(250,1000)
        # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.png')
        # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase.pdf')
        print("--- %s seconds ---" % (time.time() - start_time))