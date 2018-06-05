#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 16:58:14 2018

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
from glob import glob

f, ax = plt.subplots(1,1)

base_dir = '/media/Raid/simulations/spherical_counterflow/ANALYSIS/new_series/RUN2/'
files = glob(path.join(base_dir, '*0.10mms*.txt'))
files.sort()

for fn in files:
    d = np.loadtxt(fn);
    Tstr = path.split(fn)[1][1:4]
    T = float(Tstr)/100
    t = d[:,0]
    L = d[:,1]
    ax.plot(t, L, label = "{:.2f} K".format(T))

ax.legend(loc='best')
ax.set_xlabel("$t$ (s)")
ax.set_ylabel("total vortex length (cm)")
f.tight_layout()
f.savefig("L_t.pdf")