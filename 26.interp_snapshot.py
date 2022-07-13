#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np

import flac
import flac_interpolate as fi
import flac_gravity3 as fg
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt

#-------------------------------------------------------------------
model = sys.argv[1]
frame = int(sys.argv[2])
#model = 'cH1404'
#frame = 120
plt.rcParams["font.family"] = "Times New Roman"
path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
#path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
#savepath='/Volumes/SSD500/data/'
figpath='/home/jiching/geoflac/figure/'
#figpath='/Users/ji-chingchen/OneDrive - 國立台灣大學/年會/2022/POSTER/'
os.chdir(path+model)

fl = flac.Flac()
x, z = fl.read_mesh(frame)

phasein     = 1
tin         = 0
vis         = 1
gravity     = 0
cpp         = 0
figure_plot = 0 # 8 figure
figure_plot1= 0 # single step viscosity
figure_plot2= 0 # single step phase 
figure_plot3= 1 # 30,60,90,150 step viscosity 
figure_plot4= 0 # 30,60,90,150 step phase
# -------------------------------------------------------------------
# domain bounds
left = -300
right = 800
up = 10
down = -200
dx = 1
dz = 1

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

def interpolate_vis(frame, xtrench):
   # domain bounds in km
   fi.xmin = xtrench + left
   fi.xmax = xtrench + right
   fi.zmin = down
   fi.zmax = up

   # resolution in km
   fi.dx = dx
   fi.dz = dz

   vx, vz, visc = fi.interpolate(frame, 'visc')
   return vx, vz, visc

itrench = find_trench_index(z)
xtrench = x[itrench,0]

###############################################

# get interpolated phase either from previous run or from original data
phasefile = 'intp3-phase.%d' % frame
if phasein:
    xx, zz, ph = interpolate_phase(frame, xtrench)
    f = open(phasefile, 'w')
    f.write('%d %d\n' % xx.shape)
    flac.printing(xx, zz, ph, stream=f)
    f.close()

# get interpolated T either from previous run or from original data
tfile = 'intp3-T.%d' % frame
if tin:
    T = fl.read_temperature(frame)
    f = open(tfile, 'w')
    f.write('%d %d\n' % x.shape)
    flac.printing(x, z, T, stream=f)
    f.close()

visfile = 'intp3-visc.%d' % frame
if vis:
    vx, vz, visc = interpolate_vis(frame, xtrench)
    f = open(visfile, 'w')
    f.write('%d %d\n' % vx.shape)
    flac.printing(vx, vz, visc, stream=f)
    f.close()
gfile = 'topo-grav.%d' % frame
if gravity:
# get topography and gravity at uniform spacing
    px, topo, topomod, gravity = fg.compute_gravity2(frame)
    # convert to km and mGal before saving
    px *= 1e-3
    topo *= 1e-3
    topomod *= 1e-3
    gravity *= 1e5
    f = open(gfile, 'w')
    flac.printing(px, topo, gravity, topomod, stream=f)
    f.close()

###############
if cpp:
    cpcmd = ''' awk '{print $1,$2,$3+0}' %(phasefile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(phasefile)s.txt
    #awk '{print $1,$2,$3+0}' %(tfile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(tfile)s.txt
    #awk '{print $1,$2,$3+0}' %(gfile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(gfile)s.txt
    awk '{print $1,$2,$3+0}' %(visfile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(visfile)s.txt
    #cp %(tfile)s  %(model)s_%(tfile)s.txt
    #cp %(gfile)s  %(model)s_%(gfile)s.txt
    mv  %(model)s_%(phasefile)s.txt %(model)s_%(visfile)s.txt /home/jiching/geoflac/data/.
    #mv  %(model)s_%(phasefile)s.txt %(model)s_%(gfile)s.txt  %(model)s_%(tfile)s.txt %(model)s_%(visfile)s.txt /home/jiching/geoflac/data/.
''' % locals()
    os.system(cpcmd)

def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)
#------------------------------------------------------------------------------
if figure_plot:
    #fig, (ax)= plt.subplots(1,2,figsize=(27,5))
    fig, (ax)= plt.subplots(1,2,figsize=(24,5))
    cc = plt.cm.get_cmap('jet')
    frame_list=[30,60,90,150]
    # frame_list=[150]
    for qq,frame in enumerate(frame_list):
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)
        bwith = 4
        #--------------------- phase plotting -------------------------
        filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
        x,z,ph=np.loadtxt(filepath).T
        ax[qq,0].scatter(x,-z,c = ph,cmap = phase19,vmax=19,vmin=1,s=1.5)
        #------------------- viscosity plotting -----------------------
        filepath = savepath+model+'_intp3-visc.'+str(frame)+'.txt'
        x,z,vis=np.loadtxt(filepath).T
        ax[qq,1].scatter(x,-z,c=vis,cmap=cc,vmin=20, vmax=27,s=1.5)
        #---------------------- plot setting --------------------------
        for uu in range(len(ax[0])):
            ax[qq,uu].set_aspect('equal')
            ax[qq,uu].contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
            ax[qq,uu].spines['bottom'].set_linewidth(bwith)
            ax[qq,uu].spines['top'].set_linewidth(bwith)
            ax[qq,uu].spines['right'].set_linewidth(bwith)
            ax[qq,uu].spines['left'].set_linewidth(bwith)
            ax[qq,uu].tick_params(axis='x', labelsize=23)
            ax[qq,uu].tick_params(axis='y', labelsize=23)
            ymajor_ticks = np.linspace(200,0,num=5)
            ax[qq,uu].set_yticks(ymajor_ticks)
            xmajor_ticks = np.linspace(250,1000,num=6)
            ax[qq,uu].set_xticks(xmajor_ticks)
            ax[qq,uu].set_xlim(250,1000)
            ax[qq,uu].set_ylim(200,-30)
            xmajor_ticks = np.linspace(500,1000,num=6)
            ax[qq,uu].set_xticks(xmajor_ticks)
            ax[qq,uu].set_xlim(450,950)
            if qq != 3:
                ax[qq,uu].axes.xaxis.set_visible(False)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase&vis_all.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_phase&vis.pdf')
#--------------------------------------------------------------------------
if figure_plot1:
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    cc = plt.cm.get_cmap('jet')
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #------------------- viscosity plotting -----------------------
    filepath = savepath+model+'_intp3-visc.'+str(frame)+'.txt'
    x,z,vis=np.loadtxt(filepath).T
    ax.scatter(x,-z,c=vis,cmap=cc,vmin=20, vmax=27,s=1)
    #---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
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
    fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis.png')
    fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis.pdf')
if figure_plot2:
    fig, (ax)= plt.subplots(1,1,figsize=(12,5))
    cc = plt.cm.get_cmap('jet')
    xt,zt = fl.read_mesh(frame)
    temp = fl.read_temperature(frame)
    bwith = 3
    #--------------------- phase plotting -------------------------
    filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
    x,z,ph=np.loadtxt(filepath).T
    ax.scatter(x,-z,c=ph,cmap=phase19,vmin=1, vmax=19,s=1)
    #---------------------- plot setting --------------------------
    ax.set_aspect('equal')
    ax.contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
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
if figure_plot3:
    # fig, (ax)= plt.subplots(3,1,figsize=(13,17))
    fig, (ax)= plt.subplots(3,1,figsize=(17,16))
    cc = plt.cm.get_cmap('jet')
    bwith = 3
    frame_list=[30,60,120]
    #--------------------- viscosity plotting -------------------------
    for qq,frame in enumerate(frame_list):
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)    
        filepath = savepath+model+'_intp3-visc.'+str(frame)+'.txt'
        x,z,vis=np.loadtxt(filepath).T
        ax[qq].scatter(x,-z,c=vis,cmap=cc,vmin=20, vmax=27,s=1)
        ax[qq].scatter(x[vis<=20.3],-z[vis<=20.3],c='w',s=1.5)
    #---------------------- plot setting --------------------------
        ax[qq].set_aspect('equal')
        ax[qq].contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=26)
        ax[qq].tick_params(axis='y', labelsize=26)
        ymajor_ticks = np.linspace(200,0,num=5)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=6)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(250,1000)
        ax[qq].set_ylim(200,-30)
        # ax[0].set_title(model,fontsize=25)
        # xmajor_ticks = np.linspace(500,1000,num=6)
        # ax[qq].set_xticks(xmajor_ticks)
        # ax[qq].set_xlim(450,950)
        if qq != 2:
            ax[qq].axes.xaxis.set_visible(False)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis_all.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis_all.pdf')
if figure_plot4:
    fig, (ax)= plt.subplots(4,1,figsize=(13,16))
    cc = plt.cm.get_cmap('jet')
    bwith = 3
    frame_list=[30,60,90,150]
    #--------------------- viscosity plotting -------------------------
    for qq,frame in enumerate(frame_list):
        filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
        x,z,ph=np.loadtxt(filepath).T
        ax[qq].scatter(x,-z,c=ph,cmap=phase19,vmin=1, vmax=19,s=1)
        xt,zt = fl.read_mesh(frame)
        temp = fl.read_temperature(frame)    
    #---------------------- plot setting --------------------------
        ax[qq].set_aspect('equal')
        ax[qq].contour(xt,-zt,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
        ax[qq].spines['bottom'].set_linewidth(bwith)
        ax[qq].spines['top'].set_linewidth(bwith)
        ax[qq].spines['right'].set_linewidth(bwith)
        ax[qq].spines['left'].set_linewidth(bwith)
        ax[qq].tick_params(axis='x', labelsize=23)
        ax[qq].tick_params(axis='y', labelsize=23)
        ymajor_ticks = np.linspace(200,0,num=5)
        ax[qq].set_yticks(ymajor_ticks)
        xmajor_ticks = np.linspace(250,1000,num=6)
        ax[qq].set_xticks(xmajor_ticks)
        ax[qq].set_xlim(400,1000)
        ax[qq].set_ylim(200,-30)
        # xmajor_ticks = np.linspace(500,1000,num=6)
        # ax[qq].set_xticks(xmajor_ticks)
        # ax[qq].set_xlim(450,950)
        if qq != 3:
            ax[qq].axes.xaxis.set_visible(False)
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis_all.png')
    # fig.savefig(figpath+model+'frame_'+str(frame)+'_interp_vis_all.pdf')
