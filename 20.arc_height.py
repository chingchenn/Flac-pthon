#!/usr/bin/env python3
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
# import creat_database as cd
import matplotlib.pyplot as plt
from Main_creat_database import oceanic_slab,nodes_to_elements
model = str(sys.argv[1])
path = '/home/jiching/geoflac/'+model+'/'
#path = '/Volumes/My Book/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
rainbow = cm.get_cmap('gray_r',end)
newcolors = rainbow(np.linspace(0, 1, end))
fig, (ax)= plt.subplots(1,1,figsize=(17,12))
#============================Time Series==========================
arc_thickness = np.zeros(end)
magma = np.zeros(end)
cc=0
arc_phase  = 14
#----------------------------Main code----------------------------
for i in range(1,end):
    phase=fl.read_phase(i)
    x,z=fl.read_mesh(i)
    ele_x, ele_z = nodes_to_elements(x,z)
    height = np.zeros(nex)
    for xx in range(0,nex):
        if phase[xx,0]==arc_phase:
            for zz in range(0,nez):
                if phase[xx,zz]!=arc_phase:
                    break
            height[xx]=z[xx,0]-z[xx,zz]
    ax.plot(ele_x[:,0],height,color=newcolors[i])
    arc_thickness[i]=max(height)
print(arc_thickness)
#----------------------------------------------------------------------------
fig.savefig('/home/jiching/geoflac/'+'figure/'+model+'_arc_thickness.png')
