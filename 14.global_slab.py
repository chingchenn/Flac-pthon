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
import scipy.optimize as so
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as f2

fig = 0
figure4=0
figure3=0
path = '/home/jiching/geoflac/figure/'
input = sys.argv[1]
#-------------------read area and get grd, trace and slab info from csv--------------------
#DIR='/home/jiching/GMT/slab/'
#ff=fs.read_data('cutmodel_'+input,DIR+input)
#area,lat,lon,az=ff.T
#line_number=len(az)
DIR='/home/jiching/GMT/slab'
#ff=fs.read_data('kkk',DIR)
ff=pd.read_csv(DIR+'/'+'kkk.csv')
for kk in range(len(ff)):
    if ff.name[kk] == input:
        break
area,rlon1,rlon2,rlat1,rlat2,grd,lon,lat,az,non,leng=ff.loc[kk].tolist()
#-------------------------------call gmt to cut the trace----------------------------------
#cmd = 'cp /home/jiching/GMT/slab/%(input)s/%(grd)s .' %locals()
#os.system(cmd)
cmd = 'gmt grdtrack -E%(lon)f/%(lat)f+a%(az)f+l%(leng)f+i0.5k -G%(grd)s>table.txt' %locals()
os.system(cmd)
cmd='''
cat table.txt | awk '{if ($3!="NaN"){print $1,$2,$3}}'> %(input)s_table4.txt
''' %locals()
os.system(cmd)
cmd='cp  %(input)s_table4.txt table.txt'%locals()
os.system(cmd)
temp=np.loadtxt('table.txt')
x,y,z=temp.T
minlat=min(y);maxlat=max(y);minlon=min(x);maxlon=max(x)
print(lat,lon,az)
#-------------------------------------call gmt to plot---------------------------------------
cmd = '''
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
gmt begin %(input)s_cross_section jpg
    gmt grdcut @earth_relief_30s -R$rlon1/$rlon2/$rlat1/$rlat2 -Gcut.nc -JM15c 
    gmt grdimage cut.nc -I+a15+ne0.5 -Cmby.cpt -t20  
    gmt basemap -LjRT+c20+w200k+f+o1c/0.2c+u -F+gwhite@50 -B
    gmt coast -W0.3p
    gmt makecpt -T-600/0/50
    gmt grdcontour $grd -C -W1p+cl -A50
    gmt colorbar -DjMR+w5c/0.3c+o-2.3c/-2.8c -Bx -By+l"km"
    gmt inset begin -DjBL+w3.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
        gmt coast  -Rg -JG$clon/$clat/? -Bg -Wfaint -G67/205/128 -A5000
        echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
    gmt inset end
    awk '{if ($3!="NaN"){print $1,$2,$3}}' table.txt > table4.txt
    mv table4.txt table.txt
    head -1 table.txt | awk '{print $1,$2}'>line.txt
    tail -1 table.txt | awk '{print $1,$2}'>>line.txt
    gmt plot -W2p,black line.txt
    gmt plot -Sc0.25c -Gblack line.txt
    gmt psbasemap -R%(minlat)f/%(maxlat)f/0/350 -JX15c/-7c -BwES -Bxa+l"Latitude (degree)" -Bya+l"Depth (km)" -Xw+3c
    awk '{print $2, (-1) * $3}' table.txt | awk '($2>0){print$1,$2}' |gmt plot -W2p
lon2=`awk '(NR==2){print $1}' line.txt`
lat2=`awk '(NR==2){print $2}' line.txt`
    gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -Gdepthtomoho.xyz > %(area)s-moho.txt
    awk '{print $2, (-1)*$4}' %(area)s-moho.txt | gmt plot -W2p,72/61/139,-- 
    gmt plot -R%(minlat)f/%(maxlat)f/-6500/5000 -Bxafg1000+l"Topography (m)" -BWsne -Bya2000f1000+l"height (m)" -JX15c/4c -W2p table.txt -Yh+0c
    gmt project -C$lon/$lat -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -Gcut.nc | awk '{print $2,$4}' >table2.txt
    gmt plot -W3p table2.txt
#    rm -f line.txt table2.txt cut.nc moho.txt
gmt end
#mv  %(input)s_cross_section* ~/geoflac/figure/.
''' %locals()
#print cmd
os.system(cmd)
