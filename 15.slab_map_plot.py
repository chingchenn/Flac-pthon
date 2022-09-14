#!/usr/bin/env python3
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
fig_GMT            = 1
fig_origin         = 1
path = '/scratch2/jiching/figure/'
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
''' %locals()
os.system(cmd)
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
gmt begin %(input)s_map png,pdf
    gmt grdcut @earth_relief_30s -R$rlon1/$rlon2/$rlat1/$rlat2 -Gcut.nc -JM15c 
#    gmt grdimage cut.nc -I+a15+ne0.5 -Cmby.cpt -t20  
    gmt grdimage cut.nc -I+a15+ne0.1  -Ccolombia.cpt   
    gmt basemap -LjRT+c20+w500k+f+o1c/0.2c+u -F+gwhite@50 -B
    gmt coast -W0.3p
    gmt makecpt -T-600/0/50
    gmt grdcontour $grd -C -W1p+cl -A50
    #gmt colorbar -DjMR+w5c/0.3c+o-2.3c/-2.8c -Bx -By+l"km"
    gmt inset begin -DjBL+w2.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
        gmt coast  -Rg -JG$clon/$clat/? -Bg -Wfaint -G67/205/128 -A5000
        echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
    gmt inset end
    head -1 table.txt | awk '{print $1,$2}'>line.txt
    tail -1 table.txt | awk '{print $1,$2}'>>line.txt
    gmt plot -W2p,black line.txt
gmt end 
gmt begin %(input)s_cross_section png,pdf
    gmt psbasemap -R%(minlat)f/%(maxlat)f/0/350 -JX15c/-7c -BwES -Bxa+l"Latitude (degree)" -Bya+l"Depth (km)" -Xw+3c # make depth-latitude plot
    awk '{print $2, (-1) * $3}' table.txt | awk '($2>0){print$1,$2}' |gmt plot -W2p
lon2=`awk '(NR==2){print $1}' line.txt`
lat2=`awk '(NR==2){print $2}' line.txt`
    gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -Gdepthtomoho.xyz > %(area)s-moho.txt
    awk '{print $2, (-1)*$4}' %(area)s-moho.txt | gmt plot -W2p,72/61/139,-- 
    gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -GCAM2016Litho.nc > %(area)s-litho.txt
    awk '{print $2, $4}' %(area)s-litho.txt | gmt plot -W2p,240/128/128,-- 
    gmt legend -DjLB+w4.0c+o0.5c -F+p1p+gbeige <<- EOF
S 0.5c - 0.9c - 2p,black 1.2c  Slab (Slab 2.0)
S 0.5c - 0.9c - 2p,240/128/128 1.2c  LAB (LITHO 1.0)
S 0.5c - 0.9c - 2p,72/61/139 1.2c  Moho (CRUST 1.0)
EOF
    gmt plot -R%(minlat)f/%(maxlat)f/-5000/5000 -Bx+l"Topography (m)" -BWsne -Bya2500+l"height (m)" -JX15c/4c -W2p table.txt -Yh+0c
    gmt project -C$lon/$lat -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -Gcut.nc | awk '{print $2,$4}' >table2.txt
    gmt plot -W3p table2.txt
    rm -f table.txt line.txt cut.nc CAM2016Litho.nc depthtomoho.xyz
gmt end
mv  %(input)s_map*  %(input)s_cross_section* ~/geoflac/figure/.
''' %locals()
if fig_GMT:
    os.system(cmd)
#==================================================================================================
temp2=np.loadtxt(str(input)+'_table4.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
x,y,z = data.T
sx = x[0]
sy = y[0]
new_cord=np.zeros(len(x))
for uu in range(1,len(x)):
    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
x=new_cord[z>mindepth]
z=z[z>mindepth]
x=x[z<-10]
z=z[z<-10]
new_cord=x

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

cmd = 'rm %(grd)s table2.txt %(input)s_table4.txt  %(area)s-moho.txt %(area)s-litho.txt' %locals()
#os.system(cmd)

