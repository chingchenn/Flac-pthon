#!/usr/bin/env python3
import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
# import creat_database as cd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
# model = str(sys.argv[1])
# path = '/home/jiching/geoflac/'+model+'/'
# path = model
# print(model)
# os.chdir(path)
model = 'w0901'
path = '/Volumes/My Book/model/'+model+'/'
print(model)
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1

melt = np.zeros(end)
magma = np.zeros(end)
cc=0
for i in range(1,end):
    mm=fl.read_fmelt(i)
    chamber=fl.read_fmagma(i)
    melt[i]=np.max(mm)
    magma[i]=np.max(chamber)
    if magma[i] >= 0.01:
        cc += 1

print("-------------------")
print(cc)
print("-------------------")
print('melt=',np.max(melt))
print(fl.time[np.argmax(melt)])
print("-------------------")
print('magma=',np.max(magma))
print(fl.time[np.argmax(magma)])
print("-------------------")

#----------------------------------------------------------------------------
def read_time(start_vts,model_steps):
    timestep=[0]
    for step in range(start_vts,model_steps+1):
        timestep.append(fl.time[step])
    # timestep=np.array(timestep)
    return timestep
def get_magma(start_vts=1,model_steps=end-1):
    melt=np.zeros(end)
    magma=np.zeros(end)
    yymelt=np.zeros(end)
    yychamber=np.zeros(end)
    arc_vol=np.zeros(end)
    for i in range(1,end):
        x,z=fl.read_mesh(i)
        phase = fl.read_phase(i)
        mm=fl.read_fmelt(i)
        chamber=fl.read_fmagma(i)
        melt[i] = np.max(mm)
        magma[i] = np.max(chamber)
        arc_vol[i]=np.sum(fl.read_area(i)[phase ==14])/1e6
        yymelt[i]=(fl.read_fmelt(i)*fl.read_area(i)/1e6).sum()
        yychamber[i]=(fl.read_fmagma(i)*fl.read_area(i)/1e6).sum()
    return melt,magma,yymelt,yychamber,arc_vol
tttime=read_time(1,end-1)
melt,magma,yymelt,yychamber,arc_vol = get_magma(1,end)
lam0 = 1e-13
lam_tdep = 8e-3
delT=40
lam = lam0*(1+np.exp(lam_tdep*delT))
prod = 6e-13
dt = 2
total_magma=prod * dt - magma * np.exp(-lam * dt) 
#----------------------------------------------------------------------------
plt.plot(tttime,total_magma)
plt.savefig('/Users/ji-chingchen/Desktop/flac-testtingcode/test.png')