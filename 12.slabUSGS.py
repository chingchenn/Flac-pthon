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

fig =1
figure2=1
figure4=1
figure3=1
residual=1

path = '/home/jiching/geoflac/figure/'
input = 'g'
#-------------------read area and get grd, trace and slab info from csv--------------------
# DIR='/home/jiching/GMT/slab'
# ff=fs.read_data('kkk',DIR)
# for kk in range(len(ff)):
#     if ff[kk,0] == input:
#         break
# area,rlon1,rlon2,rlat1,rlat2,grd,lon1,lon2,lat1,lat2=ff[kk,:]
# maxlon = max(lon1,lon2)
# minlon = min(lon1,lon2)
# maxlat = max(lat1,lat2)
# minlat = min(lat1,lat2)
# #print(f2.getDistance(lat1,lon1,lat2,lon2))
# clon = (lon1+lon2)/2
# clat = (lat1+lat2)/2
# #-------------------------------call gmt to cut the trace----------------------------------
# cmd = 'cp /home/jiching/GMT/slab/depgrd/%(grd)s .' %locals()
# os.system(cmd)
# cmd = 'gmt grdtrack -E%(lon1)f/%(lat1)f/%(lon2)f/%(lat2)f+i0.5k -G%(grd)s >table.txt' %locals()
# os.system(cmd)
# #-------------------------------------call gmt to plot---------------------------------------
# cmd = '''
# grd=%(grd)s
# rlon1=%(rlon1)f
# rlon2=%(rlon2)f
# rlat1=%(rlat1)f
# rlat2=%(rlat2)f
# lon1=%(lon1)f
# lat1=%(lat1)f
# lon2=%(lon2)f
# lat2=%(lat2)f
# clon=`echo %(rlon1)f %(rlon2)f| awk '{print ($1 + $2)/2 }'`
# clat=`echo %(rlat1)f %(rlat2)f| awk '{print ($1 + $2)/2 }'`
# gmt begin %(input)s_cross_section jpg
#     gmt grdcut @earth_relief_30s -R$rlon1/$rlon2/$rlat1/$rlat2 -GcutMexico.nc -JM15c 
#     gmt grdimage cutMexico.nc -I+a15+ne0.5 -Cmby.cpt -t20 --FORMAT_GEO_MAP=dddF 
#     gmt basemap -LjRT+c20+w200k+f+o1c/0.2c+u -F+gwhite@50 -B
#     gmt coast -W0.3p
#     gmt makecpt -T-600/0/50
#     gmt grdcontour $grd -C -W1p+cl -A50
#     gmt colorbar -DjMR+w5c/0.3c+o-2.3c/-2.8c -Bx -By+l"km"
#     cat <<- EOF > line.txt
#     $lon1   $lat1
#     $lon2   $lat2
# EOF
#     gmt plot -W2p,black line.txt
#     gmt inset begin -DjBL+w3.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
#         gmt coast  -Rg -JG%(clon)f/%(clat)f/? -Bg -Wfaint -G67/205/128 -A5000
#         echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
#     gmt inset end
#     gmt plot -Sc0.25c -Gblack line.txt
#     gmt grdtrack -E$lon1/$lat1/$lon2/$lat2+i0.5k -G$grd >table.txt
#     gmt psbasemap -R%(minlat)f/%(maxlat)f/0/350 -JX15c/-7c -BwES -Bxa+l"Latitude (degree)" -Bya+l"Depth (km)" -Xw+3c
#     #awk '{print $2, (-1) * $3}' table.txt |gmt plot -W2p
#     awk '{print $2, (-1) * $3}' table.txt | awk '($2>0){print$1,$2}' |gmt plot -W2p
#     gmt plot -R%(minlat)f/%(maxlat)f/-6500/5000 -Bxafg1000+l"Topography (m)" -BWsne -Bya2000f1000+l"height (m)" -JX15c/4c -W2p table.txt -Yh+0c
#     gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -GcutMexico.nc | awk '{print $2,$4}' >table2.txt
#     gmt plot -W3p table2.txt
#     rm -f line.txt table2.txt cutMexico.nc
# gmt end
# mv  %(input)s_cross_section* ~/geoflac/figure/.
# ''' %locals()
#print cmd
# os.system(cmd)
# #------------------------------------------------------------------------------------------
path='/Users/ji-chingchen/OneDrive - 國立台灣大學/resarch/slab 2.0/slab/'
file='Gline.txt'
temp2=np.loadtxt(path+file)
# temp2 = np.loadtxt('S.txt')
data = temp2[~np.isnan(temp2).any(axis=1)]
x,y,z = data.T
sx = x[0]
sy = y[0]
new_cord=np.zeros(len(x))
for uu in range(1,len(x)):
    new_cord[uu]=f2.getDistance(y[uu], x[uu], sy, sx)
def find_Ct(A):
    C=np.zeros((len(A),len(A[0])))
    for kk in range(len(A)):
        new_array = np.delete(A,kk,axis=0)
        for qq in range(len(A[0])): 
            new_array_2 = np.delete(new_array,qq,axis=1)
            ww=np.linalg.det(new_array_2)
            C[kk][qq] = ww*(-1)**(kk+qq) 
            cc=C.T
    return cc
def find_inv(A):
    adjA = find_Ct(A)
    detA = abs(np.linalg.det(A))
    return adjA/detA
# plt.plot(new_cord,z)
# plt.savefig(path+input+'wwcross+poly.png')
xx=new_cord
zz=z
mindepth=-250
x=new_cord[z>mindepth]
z=z[z>mindepth]
N=len(x)
G = np.array([np.ones(N),x])
GT=G.T
#m1=find_inv(GT.dot(G)).dot(GT).dot(x)

# Polynomail 4
z4=np.polyfit(x,z,4)
w4=np.polyval(z4,x)
res4=sum((w4-z)**2)
sst=sum((z-np.mean(z))**2)
R4=1-(res4/sst)

# Polynomial 3
z3=np.polyfit(x,z,3)
w3=np.polyval(z3,x)
res3=sum((w3-z)**2)
R3=1-(res3/sst)

# Polynomial 2
z2=np.polyfit(x,z,2)
w2=np.polyval(z2,x)
res2=sum((w2-z)**2)
R2=1-(res2/sst)

# Polynomial 1
z1=np.polyfit(x,z,1)
w1=np.polyval(z1,x)
res1=sum((w1-z)**2)
R1=1-(res1/sst)


rr=[R1,R2,R3,R4]
nn=[1,2,3,4]
kk=['x^2','x^3','x^4']

if fig:
    fig, (ax)= plt.subplots(1,1,figsize=(10,8))
    ax.plot(x,w4,c='#4169E1',lw=3,label='quartic')
    ax.plot(x,w3,c='r',lw=2,label='cubic')
    ax.plot(x,w2,c='orange',lw=2,label='quadratic')
    ax.plot(x,w1,c='green',lw=2,label='linear')
    ax.grid()
    ax.set_aspect('equal', adjustable='box')
    # ax.set_ylim(mindepth,0)
    ax.plot(x,z,'k--',label='observation') 
    ax.plot(xx,zz,'k--')#,label='observation') 
    # ax.set_ylim(-300,0)
    # ax.set_xlim(0,600)
    ax.legend()
    # fig.savefig(path+input+'cross+poly.png')

if figure2:
    fig2, ax2 = plt.subplots(1,1,figsize=(6,8))
    ax2.scatter(nn,rr,s=50,c='r')
    ax2.set_ylim(0.7,1.1)
    ax2.set_xlim(0,5)
    ax2.grid()
    ax2.set_ylabel('R^2',fontsize=20)
    # fig2.savefig(path+input+'rsquare.png')

if figure4:
    p4=np.poly1d(z4)
    fp3=np.polyder(p4,1)
    fp2=np.polyder(p4,2)
    f3=fp3(x)
    f2=fp2(x)
    fig3, (ax3,ax4,ax5)= plt.subplots(3,1,figsize=(9,12))
    ax3.plot(x,z,c='#4169E1',lw=2)
    ax3.plot(x,w4,c='k')
    ax4.plot(x,f3,c='k')
    ax5.plot(x,f2,c='k')
    ax3.grid();ax4.grid();ax5.grid()
    ax3.set_title('quartic',fontsize=20)
    # fig3.savefig(path+input+'poly4_analyses.png')
if figure3:
    p3=np.poly1d(z3)
    fp3=np.polyder(p3,1)
    # fp2=np.polyder(p3,2)
    f3=fp3(x)
    # f2=fp2(x)
    fig3, (ax3,ax4)= plt.subplots(2,1,figsize=(9,8))
    ax3.plot(x,z,c='#4169E1',lw=2)
    ax3.plot(x,w4,c='k')
    ax4.plot(x,f3,c='k')
    # ax5.plot(x,f2,c='k')
    ax3.grid();ax4.grid()
    ax3.set_title('cubic',fontsize=20)
    
if residual:
    r4=(w4-z)
    fig5, (ax7,ax8,ax9,ax10)= plt.subplots(4,1,figsize=(9,12))
    ax7.scatter(x,r4,c='b',s=2)
    r3=w3-z
    ax8.scatter(x,r3,c='b',s=2)
    r2=w2-z
    ax9.scatter(x,r2,c='b',s=2)
    ax7.set_title('Residual',fontsize=20)
    r1=w1-z
    ax10.scatter(x,r1,c='b',s=2)

#cmd = 'rm %(grd)s table.txt' %locals()
#os.system(cmd)
