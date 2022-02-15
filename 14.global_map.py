#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:14:24 2021
@author: jiching
"""
import flac
import os,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as f2

fig = 0
input='kk'
path = '/home/jiching/geoflac/figure/'
#-------------------------------call gmt to cut the trace----------------------------------
cmd = 'cp /home/jiching/GMT/slab/depgrd/*grd .' %locals()
os.system(cmd)
#cmd = 'gmt grdtrack -E%(lon1)f/%(lat1)f/%(lon2)f/%(lat2)f+i0.5k -G%(grd)s >table.txt' %locals()
#os.system(cmd)
#-------------------------------------call gmt to plot---------------------------------------
cmd = '''

gmt begin All png
    gmt coast -JH180/12c -Rg -Bg -W0.5p -A10000
    gmt grdcut @earth_relief_20m -Rg -Gcut.nc
    gmt grdimage cut.nc -I+a15+ne0.3 -Cgray.cpt -t20 -JM15c --FORMAT_GEO_MAP=dddF
    gmt basemap -LjRT+c20+w5000k+f+o1c/0.2c+u -F+gwhite@50 -B
    gmt coast -W0.3p
    gmt grdcontour alu_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour cal_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour cam_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour car_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour cas_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour cot_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour hal_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour hel_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour him_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour hin_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour izu_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour ker_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour kur_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour mak_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour man_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour mue_slab2_dep_02.24.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour phi_slab2_dep_02.26.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour pam_slab2_dep_02.26.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour png_slab2_dep_02.26.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour puy_slab2_dep_02.26.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour ryu_slab2_dep_02.26.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour sam_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour sco_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour sol_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour sul_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour sum_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    gmt grdcontour van_slab2_dep_02.23.18.grd -Crainbow.cpt -W1p+cl -A20
    awk -F, '(NR>1){print $7, $8,$9,$10}' ../GMT/slab/kkk.csv|gmt plot -Sc0.2c -G#006400
rm cut.nc
gmt end
#mv  %(input)s_cross_section* ~/geoflac/figure/.
''' %locals()
#print cmd
os.system(cmd)
# #------------------------------------------------------------------------------------------
#temp2 = np.loadtxt('trenchamz.txt')
#DIR='/home/jiching/GMT/slab/trenchaz'
#DIR='/home/jiching/GMT/slab'
#ff=fs.read_data('trenchaz',DIR)
#for kk in range(len(ff)):
#    name,lat,lon,az=ff[kk].T
        

#data = temp2[~np.isnan(temp2).any(axis=1)]
#x,y,z = data.T
#sx = x[0]
#isy = y[0]
#new_cord=np.zeros(len(x))
#for uu in range(1,len(x)):
#    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
#z1=np.polyfit(new_cord,z,4)
#p4=np.poly1d(z1)
#w1=p4(new_cord)
#fig, (ax)= plt.subplots(1,1,figsize=(10,12))
#ax.scatter(x,z,s=2,color='#4169E1')
#ax.plot(x,w1,lw=2,color='k')
#ax.set_ylim(-300,0)
#ax.set_xlim(-100.25,-98.35)

#if fig:
#    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
#    ax.scatter(new_cord,z,s=2,color='#4169E1')
#    ax.set_aspect('equal', adjustable='box')
#    ax.set_ylim(-300,0)
#    ax.set_xlim(0,650)
#    fig.savefig(path+input+'_cross.png')
#cmd = 'rm %(grd)s table.txt' %locals()
#os.system(cmd)
