#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 11:33:06 2018

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
from glob import glob

base_dir = '/media/Raid/simulations/spherical_counterflow/ANALYSIS/new_series/RUN2/'
files = glob(path.join(base_dir, '*L_t.txt'))
files.sort()

Ts = []
vs = []
Ls = []


for fn in files:
    d = np.loadtxt(fn);
    basename = path.split(fn)[1]
    Tstr = basename[1:4]
    vstr = basename[5:-11]
    T = float(Tstr)/100
    v = float(vstr)
    L = d[:,1][-1]

    Ts.append(T)
    vs.append(v)
    Ls.append(L)

Ls = np.array(Ls)

f, ax = plt.subplots(1,1)
ax.scatter(Ts, vs, sizes = 10*Ls)
ax.set_yscale('log')

ax.set_xlabel("$T$ (K)")
ax.set_ylabel(r"$v_\mathrm{ns}$ at 5 mm (mm/s)")
f.tight_layout()
f.savefig("phase_diagram.pdf")