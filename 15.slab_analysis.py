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
fig_origin         = 0
fig_spline         = 0
fig_poly           = 0
fig_Rsquare        = 0
fig_quartic        = 0
fig_cubic          = 0
fig_residual       = 0
fig_spline_quartic = 0
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
#    gmt colorbar -DjMR+w5c/0.3c+o-2.3c/-2.8c -Bx -By+l"km"
    gmt inset begin -DjBL+w3.2c+o0.3c/0.3c -F+gwhite+p1p+c0.1c
        gmt coast  -Rg -JG$clon/$clat/? -Bg -Wfaint -G67/205/128 -A5000
        echo $rlon1 $rlat1 $rlon2 $rlat2 | gmt plot -Sr+s -W1p,blue
    gmt inset end
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
    gmt project -C$lon/$lat -E$lon2/$lat2 -Q -G0.1 | gmt grdtrack -GCAM2016Litho.nc > %(area)s-litho.txt
    awk '{print $2, $4}' %(area)s-litho.txt | gmt plot -W2p,240/128/128,-- 
    gmt legend -DjLB+w5.5c+o0.5c -F+p1p+gbeige <<- EOF
S 0.5c - 0.9c - 2p,black 1.2c  Slab 2.0
S 0.5c - 0.9c - 2p,240/128/128 1.2c  LAB (LITHO 1.0)
S 0.5c - 0.9c - 2p,72/61/139 1.2c  Moho (CRUST 2.0)
EOF
    gmt plot -R%(minlat)f/%(maxlat)f/-6500/5000 -Bxafg1000+l"Topography (m)" -BWsne -Bya2000f1000+l"height (m)" -JX15c/4c -W2p table.txt -Yh+0c
    gmt project -C$lon/$lat -E$lon2/$lat2 -G0.1 -Q | gmt grdtrack -Gcut.nc | awk '{print $2,$4}' >table2.txt
    gmt plot -W3p table2.txt
    rm -f line.txt table2.txt cut.nc %(area)s-moho.txt table.txt %(area)s-litho.txt CAM2016Litho.nc depthtomoho.xyz
gmt end
mv  %(input)s_cross_section* ~/geoflac/figure/.
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
#==================================================BSpline==========================================
kk=3
ss=0.001
mean_list=np.zeros(len(x))
median_list=np.zeros(len(x))
tck = interpolate.splrep(x,z,k=kk,s=ss)
zz0=interpolate.splev(x,tck,der=0)    
zz1=interpolate.splev(x,tck,der=1)    
zz2=interpolate.splev(x,tck,der=2) 
yders = interpolate.spalde(x, tck)
#print('k=',kk,'s=',ss,'mean=',np.mean(zz0-z),'median=',np.median(zz0-z))
if fig_spline:
    fig0,(q1)= plt.subplots(1,1,figsize=(10,8))
    q1.plot(x,zz0,c='r')
    q1.plot(new_cord,z,'k--')
    q1.set_aspect('equal', adjustable='box')
#    q2.plot(x,zz1,c='r')
#    q3.plot(x,zz2,c='r')
    q1.set_ylim(mindepth,0)
    q1.grid()
#    q2.grid();q3.grid()
    q1.tick_params(axis='x', labelsize=16)
    q1.tick_params(axis='y', labelsize=16)
#    q2.tick_params(axis='x', labelsize=16)
#    q2.tick_params(axis='y', labelsize=16)
#    q3.tick_params(axis='y', labelsize=16)
#    q3.tick_params(axis='x', labelsize=16)
    fig0.savefig(path+str(input)+'_slab_spline.png')

#==============================================Polynomail============================================
## Polynomail 4
z4=np.polyfit(x,z,4)
w4=np.polyval(z4,x)
res4=sum((w4-z)**2)
sst=sum((z-np.mean(z))**2)
R4=1-(res4/sst)

## Polynomial 3
z3=np.polyfit(x,z,3)
w3=np.polyval(z3,x)
res3=sum((w3-z)**2)
R3=1-(res3/sst)

## Polynomial 2
z2=np.polyfit(x,z,2)
w2=np.polyval(z2,x)
res2=sum((w2-z)**2)
R2=1-(res2/sst)

## Polynomial 1
z1=np.polyfit(x,z,1)
w1=np.polyval(z1,x)
res1=sum((w1-z)**2)
R1=1-(res1/sst)

RR=[R1,R2,R3,R4]
nn=[1,2,3,4]
#---------------------------------data result without fitting-----------------------------------
m=[]; m2=[]
for kk in range(1,len(x)):
    cx1=x[kk-1];cx2=x[kk]
    cz1=z[kk-1];cz2=z[kk]
    if (cx2-cx1) != 0:
      m.append((cz2-cz1)/(cx2-cx1))
qq = x[1:]
for ww in range(1,len(qq)):
    cx1=qq[ww-1];cx2=qq[ww]
    cz1=m[ww-1];cz2=m[ww]
    if (cx2-cx1) != 0:
        m2.append((cz2-cz1)/(cx2-cx1))
qq2=qq[1:ww+1]
#-------------------------defint flat slab or not by Bspline result ---------------------------
ff1=[];therd=0.2
#for rr,oo in enumerate(zz1): # first derivative
#    if (oo+therd)>0:
#        ff1.append(x[rr])
#    cc = oo
mm=-1;ff2=[]
for pp,uu in enumerate(zz2): # second derivative
    if mm*uu<0:
        ff2.append(x[pp])
    mm = uu
for ee in range(1,len(ff2)):
    if  zz2[(x>ff2[ee-1])* (x<ff2[ee])].all():
        cond1=ff2
        break
for anyone in (zz1[(x>cond1[0])* (x<cond1[1])]):
    if (anyone + therd)>0:
        ff1.append(anyone)
if len(ff2)>1 and (ff2[1]-ff2[0])>90 :
    print(str(input)+' has the possibility to be a flat slab')
    if len(ff1)>10:
        print(str(input)+' is flat slab')
else: print(str(input)+' is not a flat slab')
print(ff2)
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
if fig_spline:
    fig0,(q1)= plt.subplots(1,1,figsize=(10,8))
    q1.plot(x,zz0,c='r')
    q1.plot(new_cord,z,'k--')
    q1.set_aspect('equal', adjustable='box')
    q1.set_ylim(mindepth,0)
    q1.grid()
    q1.tick_params(axis='x', labelsize=16)
    q1.tick_params(axis='y', labelsize=16)
    fig0.savefig(path+str(input)+'_slab_spline.png')
if fig_poly:
    fig, (ax)= plt.subplots(1,1,figsize=(10,8))
    ax.plot(x,w4,c='#4169E1',lw=3,label='quartic')
    ax.plot(x,w3,c='r',lw=2,label='cubic')
    ax.plot(x,w2,c='orange',lw=2,label='quadratic')
    ax.plot(x,w1,c='green',lw=2,label='linear')
    ax.grid()
    ax.set_aspect('equal', adjustable='box')
    ax.plot(x,z,'k--',label='observation') 
    # ax.set_ylim(-300,0)
    # ax.set_xlim(0,600)
    ax.legend()
    ax.set_ylim(mindepth,0)
    ax.plot(x,z,color='#4169E1') 
    ax.set_ylim(-150,0)
    ax.set_xlim(0,400)
    fig.savefig(path+input+'_slab_poly.png')
if fig_Rsquare:
    fig2, ax2 = plt.subplots(1,1,figsize=(6,8))
    ax2.scatter(nn,RR,s=50,c='r')
    ax2.set_ylim(0.7,1.1)
    ax2.set_xlim(0,5)
    ax2.grid()
    ax2.set_ylabel('R^2',fontsize=20)
    fig2.savefig(path+input+'slab_rsquare.png')
if fig_quartic:
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
    fig3.savefig(path+input+'slab_poly4_analyses.png')
if fig_cubic:
    p3=np.poly1d(z3)
    fp3=np.polyder(p3,1)
    f3=fp3(x)
    fig4, (ax3,ax4)= plt.subplots(2,1,figsize=(9,8))
    ax3.plot(x,z,c='#4169E1',lw=2)
    ax3.plot(x,w4,c='k')
    ax4.plot(x,f3,c='k')
    ax3.grid();ax4.grid()
    ax3.set_title('cubic',fontsize=20)
    fig4.savefig(path+input+'slab_poly3_analyses.png')
    
if fig_residual:
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
    fig5.savefig(path+input+'slab_poly_residual.png')
if fig_spline_quartic:
    p4=np.poly1d(z4)
    fp3=np.polyder(p4,1)
    fp2=np.polyder(p4,2)
    f3=fp3(x)
    f2=fp2(x)
    fig6,(ax11,ax12,ax13)=plt.subplots(3,1,figsize=(10,15))
    ax11.plot(x,z,c='g',lw=5)
    ax11.plot(x,z,c='#4169E1',lw=5,label='data')
    ax12.plot(qq,m,c='#4169E1',lw=5,label='data')
    ax13.plot(qq2,m2,c='#4169E1',lw=5,label='data')
    ax11.plot(x,w4,c='k',label='quartic')
    ax12.plot(x,f3,c='k',label='quartic')
    ax13.plot(x,f2,c='k',label='quartic')
    ax11.grid();ax12.grid();ax13.grid()
    ax11.set_title('quartic, spline and data',fontsize=20)
    ax11.plot(x,zz0,c='r',label='Bspline')
    ax12.plot(x,zz1,c='r',label='Bspline')
    ax13.plot(x,zz2,c='r',label='Bspline')
    ax11.set_ylim(mindepth,0)
    ax11.tick_params(axis='x', labelsize=16)
    ax11.tick_params(axis='y', labelsize=16)
    ax12.tick_params(axis='x', labelsize=16)
    ax12.tick_params(axis='y', labelsize=16)
    ax13.tick_params(axis='y', labelsize=16)
    ax13.tick_params(axis='x', labelsize=16)
    fig6.savefig(path+input+'slab_ploy+spline.png')

cmd = 'rm %(grd)s %(input)s_table4.txt' %locals()
os.system(cmd)

