#!/usr/bin/env python

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

path='/home/jiching/geoflac/'
#path = '/scratch2/jiching/22winter/'
#path = '/scratch2/jiching/03model/'
#path = 'F:/model/'
# path = 'D:/model/'
path = '/Volumes/SSD500/model/'
savepath='/home/jiching/geoflac/data/'
savepath='/Volumes/SSD500/data/'
figpath='/home/jiching/geoflac/figure/'
os.chdir(path+model)

fl = flac.Flac()
x, z = fl.read_mesh(frame)

phasein = 1
tin     = 1
ph2grd  = 0
t2grd   = 0
gmtplot = 0
#-------------------------------------------------------------------
# domain bounds
left = -300
right = 800
up = 10
down = -200
dx = 1.5
dz = 0.6

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
cpcmd = '''
awk '{print $1,$2,$3+0}' %(phasefile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(phasefile)s.txt
awk '{print $1,$2,$3+0}' %(tfile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(tfile)s.txt
awk '{print $1,$2,$3+0}' %(gfile)s | awk '{if ($3>0) print $1,$2,$3}' > %(model)s_%(gfile)s.txt
#cp %(tfile)s  %(model)s_%(tfile)s.txt
#cp %(gfile)s  %(model)s_%(gfile)s.txt
mv  %(model)s_%(phasefile)s.txt %(model)s_%(gfile)s.txt  %(model)s_%(tfile)s.txt /home/jiching/geoflac/data/.
''' % locals()
os.system(cpcmd)

#------------------------------------------------------------------------------
colors = ["#93CCB1","#550A35","#2554C7","#008B8B","#4CC552",
          "#2E8B57","#524B52","#D14309","#ed45a7","#FF8C00",
          "#FF8C00","#455E45","#F9DB24","#c98f49","#525252",
          "#F67280","#00FF00","#FFFF00","#7158FF"]
phase19= matplotlib.colors.ListedColormap(colors)

fig, (ax)= plt.subplots(1,1,figsize=(20,6))
filepath = savepath+model+'_intp3-phase.'+str(frame)+'.txt'
x,z,ph=np.loadtxt(filepath).T
ax.scatter(x,z,c = ph,cmap = phase19,vmax=19,vmin=1,s=5)
ax.set_aspect('equal')
def nodes_to_elements(xmesh,zmesh):
    ele_x = (xmesh[:fl.nx-1,:fl.nz-1] + xmesh[1:,:fl.nz-1] + xmesh[1:,1:] + xmesh[:fl.nx-1,1:]) / 4.
    ele_z = (zmesh[:fl.nx-1,:fl.nz-1] + zmesh[1:,:fl.nz-1] + zmesh[1:,1:] + zmesh[:fl.nx-1,1:]) / 4.
    return ele_x, ele_z
x,z = fl.read_mesh(frame)
temp = fl.read_temperature(frame)
ax.contour(x,z,temp,cmap='rainbow',levels =[0,200,400,600,800,1000,1200],linewidths=3)
ax.set_xlim(200,1000)
ax.set_ylim(-200,10)

#--------------------------------------------------------------------------



