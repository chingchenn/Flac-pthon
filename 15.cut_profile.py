#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:14:24 2021
@author: JiChing Chen
"""
import flac
import os,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import scipy.optimize as so
from scipy import interpolate
import function_savedata as fs
import function_for_flac as f2
import matplotlib.pyplot as plt

# ===================================initial set up ======================================
cut_profile          = 1
fig_origin         = 0
path = '/home/jiching/geoflac/figure/'
input = sys.argv[1]
#-------------------read area and get grd, trace and slab info from csv--------------------
DIR='/home/jiching/GMT/slab'
ff=pd.read_csv(DIR+'/'+'kkk.csv')
for kk in range(len(ff)):
    if ff.name[kk] == input:
        break
area,rlon1,rlon2,rlat1,rlat2,grd,lon,lat,az,mindepth,leng=ff.loc[kk].tolist()
#-------------------------------call gmt to cut the trace----------------------------------
cmd = '''
cp /home/jiching/GMT/slab/depgrd/%(grd)s .
cp /home/jiching/GMT/slab/global_model/CAM2016Litho.nc .
cp /home/jiching/GMT/slab/global_model/depthtomoho.xyz .
cp /home/jiching/observation/gravity/grav_31.1.nc .
''' %locals()
os.system(cmd)
cmd = 'gmt grdtrack -E%(lon)f/%(lat)f+a%(az)f+l%(leng)f+i0.5k -G%(grd)s>slab_all.txt' %locals()
os.system(cmd)
cmd='''
cat slab_all.txt | awk '{if ($3!="NaN"){print $1,$2,$3}}'> %(input)s_slab.txt
''' %locals()
os.system(cmd)
cmd='cp  %(input)s_slab.txt slab.txt'%locals()
os.system(cmd)
temp=np.loadtxt('slab.txt')
x,y,z=temp.T
minlat=min(y);maxlat=max(y);minlon=min(x);maxlon=max(x)
#-------------------------------------call gmt to plot---------------------------------------
cmd = '''
gmt set FONT_ANNOT_PRIMARY          8p,4,black,Times-Roman \
        FONT_TITLE                  30p,4,black,Times-Roman

grd=%(grd)s
rlon1=%(rlon1)f
rlon2=%(rlon2)f
rlat1=%(rlat1)f
rlat2=%(rlat2)f
lon=%(lon)f
lat=%(lat)f
az=%(az)f
len=%(leng)f
clon=`echo %(rlon1)f %(rlon2)f| awk '{print ($1 + $2)/2 }'`
clat=`echo %(rlat1)f %(rlat2)f| awk '{print ($1 + $2)/2 }'`

gmt grdcut @earth_relief_30s -R$rlon1/$rlon2/$rlat1/$rlat2 -Gcut.nc -JM15c 
head -1 slab.txt | awk '{print $1,$2}'>line.txt
tail -1 slab.txt | awk '{print $1,$2}'>>line.txt
lon2=`awk '(NR==2){print $1}' line.txt`
lat2=`awk '(NR==2){print $2}' line.txt`
gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -Gdepthtomoho.xyz > %(area)s-moho.txt
gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -GCAM2016Litho.nc > %(area)s-litho.txt
gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -Ggrav_31.1.nc > %(area)s-grav.txt 
gmt project -C$lon/$lat -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -Gcut.nc | awk '{print $2,$4}' >%(area)s-topo.txt 

''' %locals()
if cut_profile   :
    os.system(cmd)
#==================================================================================================
temp2=np.loadtxt(str(input)+'-topo.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
#x,y,z = data.T
#sx = x[0]
#sy = y[0]
#new_cord=np.zeros(len(x))
#for uu in range(1,len(x)):
#    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
#x=new_cord[z>mindepth]
#z=z[z>mindepth]
#x=x[z<-10]
#z=z[z<-10]
#new_cord=x
#-------------------------------------call pyhton to plot---------------------------------------

if fig_origin:
    ff, (aa) = plt.subplots(1,1,figsize=(6,3))
    aa.plot(new_cord,z,'k--')
    aa.set_aspect('equal', adjustable='box')
    aa.tick_params(axis='x', labelsize=16)
    aa.tick_params(axis='y', labelsize=16)
    aa.set_ylim(-200,0)
    aa.set_xlim(0,800)
    ff.savefig(path+str(input)+'_slab.png')
