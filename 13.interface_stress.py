# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:54:38 2021

@author: jiching
"""
import flac
import math
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import function_for_flac as f2

model='s0917'
path = 'D:/model/'+model+'/'
os.chdir(path)
fl = flac.Flac();end = fl.nrec
nex = fl.nx - 1;nez = fl.nz - 1
flame = 255
strain = fl.read_strain(flame)
sxx = fl.read_sxx(flame)
strain_rate = fl.read_srII(flame)
x,y = fl.read_mesh(flame)
xx = fl.nx
zz = fl.nz
element_x = 1/4*(x[:xx-1,:zz-1]+x[1:xx,1:zz]+x[1:xx,:zz-1]+x[:xx-1,1:zz])
element_y = 1/4*(y[:xx-1,:zz-1]+y[1:xx,1:zz]+y[1:xx,:zz-1]+y[:xx-1,1:zz])
    
for kk in range(nez):
    maxsr = np.amax(strain_rate[:,kk])
    maxsr_id=np.argmin(maxsr)
    print(maxsr)
    # maxsr_id = strain_rate[:,kk].index(maxsr)
plt.scatter(element_x,element_y,c=strain_rate,s=2.5)


