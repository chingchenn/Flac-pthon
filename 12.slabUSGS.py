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
import matplotlib.pyplot as plt
import function_savedata as fs
import function_for_flac as f2

fig = 1
path = '/home/jiching/geoflac/figure/'
input = 'g'
#-------------------read area and get grd, trace and slab info from csv--------------------
DIR='/home/jiching/GMT/slab'
ff=fs.read_data('kkk',DIR)
for kk in range(len(ff)):
    if ff[kk,0] == input:
        break
area,rlon1,rlon2,rlat1,rlat2,grd,lon1,lon2,lat1,lat2=ff[kk,:]
maxlon = max(lon1,lon2)
minlon = min(lon1,lon2)
maxlat = max(lat1,lat2)
minlat = min(lat1,lat2)
print(f2.getDistance(lat1,lon1,lat2,lon2))
clon = (lon1+lon2)/2
clat = (lat1+lat2)/2
#-------------------------------call gmt to cut the trace----------------------------------
cmd = 'cp /home/jiching/GMT/slab/depgrd/%(grd)s .' %locals()
os.system(cmd)
cmd = 'gmt grdtrack -E%(lon1)f/%(lat1)f/%(lon2)f/%(lat2)f+i0.5k -G%(grd)s >table.txt' %locals()
os.system(cmd)
#-------------------------------------call gmt to plot---------------------------------------
cmd = '''
grd=%(grd)s
rlon1=%(rlon1)f
rlon2=%(rlon2)f
rlat1=%(rlat1)f
rlat2=%(rlat2)f
lon1=%(lon1)f
lat1=%(lat1)f
lon2=%(lon2)f
lat2=%(lat2)f
clon=`echo %(rlon1)f %(rlon2)f| awk '{print ($1 + $2)/2 }'`
clat=`echo %(rlat1)f %(rlat2)f| awk '{print ($1 + $2)/2 }'`
gmt begin %(input)s_cross_section jpg
    gmt coast -R$rlon1/$rlon2/$rlat1/$rlat2 -B -G67/205/108 -W1p -t20 -JM15c
    gmt grdcut @earth_relief_30s -R$rlon1/$rlat1/$rlon2/$rlat2 -GcutMexico.nc
    gmt grdimage cutMexico.nc -I+a15+ne0.5 -Cmby.cpt -t20 --FORMAT_GEO_MAP=dddF
    gmt basemap -LjRT+c20+w200k+f+o1c/0.2c+u -F+gwhite@50
    gmt makecpt -T-600/0/50
    gmt grdcontour $grd -C -W1p+cl -A50
    gmt colorbar -DjMR+w5c/0.3c+o0.3c -Bx -By+l"km"
    cat <<- EOF > line.txt
    $lon1   $lat1
    $lon2   $lat2
EOF
    gmt plot -W2p,black line.txt
    gmt inset begin -DjBL+w3.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
        gmt coast  -Rg -JG%(clon)f/%(clat)f/? -Bg -Wfaint -G67/205/128 -A5000
        echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
    gmt inset end
    gmt plot -Sc0.25c -Gblack line.txt
    gmt grdtrack -E$lon1/$lat1/$lon2/$lat2+i0.5k -G$grd >table.txt
    gmt psbasemap -R%(minlat)f/%(maxlat)f/0/350 -JX15c/-7c -BwES -Bxa+l"Longitude (degree)" -Bya+l"Depth (km)" -Xw+3c
    awk '{print $2, (-1) * $3}' table.txt |gmt plot -W2p
    gmt plot -R%(minlat)f/%(maxlat)f/-6000/5000 -Bxafg1000+l"Topography (m)" -BWsne -JX15c/4c -W2p table.txt -Yh+0c
    gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -GcutMexico.nc | awk '{print $2,$4}' >table2.txt
    gmt plot -W3p table2.txt
    rm -f line.txt table2.txt cutMexico.nc
gmt end
mv  %(input)s_cross_section* ~/geoflac/figure/.
''' %locals()
#print cmd
os.system(cmd)
# #------------------------------------------------------------------------------------------
temp2 = np.loadtxt('table.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
x,y,z = data.T
sx = x[0]
sy = y[0]
new_cord=np.zeros(len(x))
for uu in range(1,len(x)):
    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
z1=np.polyfit(new_cord,z,4)
p4=np.poly1d(z1)
w1=p4(new_cord)
#fig, (ax)= plt.subplots(1,1,figsize=(10,12))
#ax.scatter(x,z,s=2,color='#4169E1')
#ax.plot(x,w1,lw=2,color='k')
#ax.set_ylim(-300,0)
#ax.set_xlim(-100.25,-98.35)

if fig:
    fig, (ax)= plt.subplots(1,1,figsize=(10,12))
    ax.scatter(new_cord,z,s=2,color='#4169E1')
    ax.set_aspect('equal', adjustable='box')
    ax.set_ylim(-300,0)
    ax.set_xlim(0,450)
    fig.savefig(path+input+'_cross.png')
cmd = 'rm %(grd)s table.txt' %locals()
os.system(cmd)
