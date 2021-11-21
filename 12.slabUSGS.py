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

fig = 1
path = '/home/jiching/geoflac/figure/'
input = 'b'
#-------------------read area and get grd, trace and slab info from csv--------------------
DIR='D:/GMT/slab'
ff=fs.read_data('kkk',DIR)
for kk in range(len(ff)):
    if ff[kk,0] == input:
        break
area,rlon1,rlon2,rlat1,rlat2,grd,lon1,lon2,lat1,lat2=ff[kk,:]


#-------------------------------call gmt to cut the trace----------------------------------
cmd = 'cp /home/jiching/GMT/slab/%(grd)s .'
os.system(cmd)
cmd = 'gmt grdtrack -E%(lon1)f/%(lat1)f/%(lon2)f/$(lat2)f+i0.5k -G%(grd)s >table.txt'
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
gmt begin Mexico_cross_section jpg
        gmt coast -R$rlon1/$rlon2/$rlat1/$rlat2 -B -G67/205/108 -W1p -t20 -JM15c
        gmt grdcut @earth_relief_30s -R$rlon1/$rlat1/$rlon2/$rlat2 -GcutMexico.nc
        gmt grdimage cutMexico.nc -I+a15+ne0.5 -t20 --FORMAT_GEO_MAP=dddF
        gmt basemap -LjRT+c20+w200k+f+o1c/0.2c+u -F+gwhite@50
        gmt makecpt -T-600/0/50
        gmt grdcontour $grd -C -W1p+cl -A50
        gmt colorbar -DjMR+w5c/0.3c+o0.3c -Bx -By+l"km"
        gmt inset begin -DjBL+w3.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
            gmt coast  -Rg -JG$clon/$clat/? -Bg -Wfaint -G67/205/128 -A5000
            echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
        gmt inset end
        cat <<- EOF > line.txt
        $lon1   $lat1
        $lon2   $lat2
        EOF
        gmt plot -RcutMexico.nc -W2p,black line.txt
        gmt plot -Sc0.25c -Gblack line.txt

        gmt grdtrack -E$lon1/$lat1/$lon2/$lat2+i0.5k -G$grd >table.txt
        gmt psbasemap -R$lat1/$lat2/0/350 -JX15c/-7c -BwES -Bxa+l"Longitude (degree)" -Bya+l"Depth (km)" -Xw+3c
        awk '{print $2, (-1) * $3}' table.txt |gmt plot -W2p
        gmt plot -R$lat1/$lat2/-6000/5000 -Bxafg1000+l"Topography (m)" -BWsne -JX15c/4c -W2p table.txt -Yh+0c
        gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -GcutMexico.nc | awk '{print $2,$4}' >table2.txt
        gmt plot -W3p table2.txt
        rm -f line.txt table.txt table2.txt cutMexico.nc
gmt end
'''
os.system(cmd)
# #------------------------------------------------------------------------------------------
temp2 = np.loadtxt('table.txt')
# x,y,z = temp2.T
# z1=np.polyfit(x,z,4)
# p4=np.poly1d(z1)
# w1=p4(x)
# #p3=np.polyder(p4,1)
# if fig :
#     fig, (ax)= plt.subplots(1,1,figsize=(10,12))
#     ax.scatter(x,z,s=2,color='#4169E1')
#     ax.plot(x,w1,ls=2,color='k')
#     fig.savefig(path+'slab.png')
# cmd = 'rm cam_slab2_dep_02.24.18.grd table.txt'
# os.system(cmd)
