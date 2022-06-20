# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:48:46 2022

@author: grace
"""

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

dir = 'D:/GMT/GMT/seismicity/'
file = 'catalogo_3_9.csv'
df = pd.read_csv(dir+file)
# mask1 = (df.Longitud>-101)
# mask2 = (df.Longitud<-97)
# mask3 = (df.Latitud<19.5)
# mask4 = (df.Latitud>18.5)
mask1 = (df.Longitud>-103)
mask2 = (df.Longitud<-95)
mask3 = (df.Latitud>19.5)
mask4 = (df.Latitud<25.5)
df = df[(mask1)&(mask2)&(mask3)&(mask4)]
import matplotlib.pyplot as plt
mm = np.array(df.Profundidad)[np.array(df.Profundidad)!='menos de 1']
kk = np.ones(len(mm))
for rr,item in enumerate(mm):
    kk[rr] = float(item)
    kk[rr] = int(kk[rr])
fig, (ax) = plt.subplots(1,1,figsize=(10,6))
bins = np.linspace(0,100,21)
n, bins, patches=ax.hist(kk, rwidth=0.85, color= '#708090',bins=bins)
ax.grid()
ax.set_xlabel("Depth (km)",fontsize = 20)
ax.set_ylabel("Number of event",fontsize = 20)
# ax.set_title("Histogram")
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.set_xlim(0,100)
# plt.show()