#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 17:36:31 2018

@author: emil
"""

import numpy as np
import os.path as path
import scipy.stats as stats
import matplotlib.pyplot as plt
from glob import glob

plt.close('all')

data_dir = '../data_spherical'
files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()
f, ax = plt.subplots(1,1)
for fn in files:
    fileout = fn.replace('.dat', '_rd.png')
    if path.isfile(fileout):
        continue;

    d = np.loadtxt(fn)
    print(fn)

    dls = np.empty((0,2))

    vidx = 0
    while True:
        ix = d[:,0] == vidx
        if not any(ix):
            break
        r = d[ix, 1:4]

        dr = np.diff(r, n=0)[1:]
        rm = 0.5*(r[1:, :] + r[:-1, :])



        mp = np.sqrt(np.sum(rm**2, axis=1))
        dl = np.sqrt(np.sum(dr**2, axis=1))

        row = np.column_stack((mp, dl))

        dls = np.row_stack((dls, row))

        vidx += 1


    rmin = 0
    rmax = 1
    rs = np.linspace(rmin, rmax, 200)

    rs = np.linspace(rmin, rmax, 2000)
    w = 0.01
    ls = []
    for r in rs:
        ix = np.logical_and(dls[:,0] > r,
                            dls[:,0] < r + w)
        ls.append(np.sum(dls[ix,1]))


    ax.clear()
    ax.plot(rs, ls)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 200)
    f.savefig(fileout)