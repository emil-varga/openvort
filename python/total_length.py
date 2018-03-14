#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:34:38 2018

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os.path as path

data_dir = '../data_spherical'

files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

lengths = []
files = files[:]

for fn in files:
    print(fn)
    d = np.loadtxt(fn)
    vidx = 0
    l = 0
    while True:
        ix = d[:,0] == vidx;
        if not any(ix):
            break

        rs = d[ix, 1:4]
        dr = np.diff(rs, axis=0)
        l += np.sum(np.sqrt(dr**2))
        vidx += 1
    lengths.append(l)

f, ax = plt.subplots(1,1)
ax.plot(lengths)