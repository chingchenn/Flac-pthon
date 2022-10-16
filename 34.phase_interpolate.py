#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:15:55 2022

@author: chingchen
"""

#-------------------------------------------------------------------
import sys, os
import numpy as np

import flac
import flac_interpolate as fi
import flac_gravity3 as fg
import matplotlib
import matplotlib.pyplot as plt

###############################################

model = sys.argv[1]
# frame = int(sys.argv[2])
# model = 'Nazca_a0624'
# frame = 10

def make_phase_interploate(model,frame):
    path='/home/jiching/geoflac/'
    path = '/Users/chingchen/Desktop/model/'


    os.chdir(path+model)

    fl = flac.Flac()
    savepath='/scratch2/jiching/data/'
    savepath = '/Users/chingchen/Desktop/data/'
    figpath='/scratch2/jiching/figure/'
    figpath = '/Users/chingchen/Desktop/figure/'
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
    # domain bounds
    left = -200
    right = 1000
    up = 10
    down = -200
    dx = 1.5
    dz = 0.6
    x, z = fl.read_mesh(frame)
    
    itrench = find_trench_index(z)
    xtrench = x[itrench,0]
    
    # get interpolated phase either from previous run or from original data
    phasefile = 'intp3-phase.%d' % frame
    xx, zz, ph = interpolate_phase(frame, xtrench)
    f = open(phasefile, 'w')
    f.write('%d %d\n' % xx.shape)
    flac.printing(xx, zz, ph, stream=f)
    f.close()
    
    # get topography and gravity at uniform spacing
    px, topo, topomod, gravity = fg.compute_gravity2(frame)
    # convert to km and mGal before saving
    px *= 1e-3
    topo *= 1e-3
    topomod *= 1e-3
    gravity *= 1e5
    gfile = 'topo-grav.%d' % frame
    f = open(gfile, 'w')
    flac.printing(px, topo, gravity, topomod, stream=f)
    f.close()
    
    
    ###############
    model = os.path.split(os.getcwd())[-1]
    phgrd = model+'_phase_%d.grd' % frame
    
    xmin = xtrench + left
    xmax = xtrench + right
    zmin = down
    zmax = up
    aspect_ratio = float(up - down) / (right - left)
    width = 6.5
    height = width * aspect_ratio
    shiftz = height + 0.3
    
    # height of gravity plot
    height2 = 1.0
    
    # gravity grid
    gravgridsize = 50
    gmin = int(gravity.min() / gravgridsize - 1) * gravgridsize
    gmax = int(gravity.max() / gravgridsize + 1) * gravgridsize
    gravann = max(abs(gmin), abs(gmax))
    # topography grid
    topogridsize = 2
    tpmin = int(topo.min() / topogridsize - 1) * topogridsize
    tpmax = int(topo.max() / topogridsize + 1) * topogridsize
    topoann = max(abs(tpmin), abs(tpmax))
    # interval of temperature contours
    cint = 200
    
    if not os.path.exists(phgrd):
        cmd = '''tail -n +2 %(phasefile)s | gmt xyz2grd -G%(phgrd)s -I%(dx)f/%(dz)f -R%(xmin)f/%(xmax)f/%(zmin)f/%(zmax)f
        mv  %(phgrd)s %(savepath)s
        ''' % locals()
        #print cmd
        os.system(cmd)


frame_list=[26,51,76,101,126,150]
for frame in frame_list:
    make_phase_interploate(model,frame)
    print('---------------end of '+str(frame)+'---------------')
<<<<<<< Updated upstream


=======
    
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
