#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:14 2021

@author: ji-chingchen
"""

import math
import flac
import os,sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

model = 'Ref'
#path = '/home/jiching/geoflac/'+model+'/'
path = '/Users/ji-chingchen/Desktop/model/'+model+'/'
#path = model
print(model)
os.chdir(path)
# fl = flac.Flac()
fl = flac.FlacFromVTK();end = fl.nrec
total_xx = np.zeros(end)
for kk in range(1,end):
    xx = fl.read_sxx(kk)
    bxx = xx[0,:]
    for yy in range(len(bxx)):
        total_xx[kk] += bxx[yy]
        # print(total_xx[kk], bxx[yy])
    
plt.plot(fl.time,total_xx)