#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 21:08:45 2018

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
from glob import glob

plt.close('all')

data_dir = '../data_oneloop'
files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

vs = []
rs = []
vs_teor = []
alphap = 1.766e-2

a = 1e-8
kappa = 9.97e-4;

for fn in files:
    d = np.loadtxt(fn)
    
    x = d[:,1]
    y = d[:,2]
    z = d[:,3]
    
    vx = d[:,4]
    vy = d[:,5]
    vz = d[:,6]
    
    x0 = x.mean()
    y0 = y.mean()
    z0 = z.mean()
    
    r_mean = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2).mean()
#    vx_mean = np.sqrt(vx**2 + vy**2 + vz**2).mean()
    vx_mean = vx.mean()
    
    v_teor = kappa/4/np.pi/r_mean*(np.log(8*r_mean/a) - 0.5)*(1-alphap)
    
    rs.append(r_mean)
    vs.append(vx_mean)
    vs_teor.append(v_teor)

f, ax = plt.subplots(1,1)
ax.plot(rs, vs, '-o')
ax.plot(rs, vs_teor, '--s')

vs = np.array(vs)
vs_teor = np.array(vs_teor)

f, ax = plt.subplots(1,1)
ax.plot(rs, 2*np.abs(vs_teor - vs)/(vs_teor + vs), '-x')