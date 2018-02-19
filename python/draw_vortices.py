#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 14:42:05 2017

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

def draw_vortices(fn, plot_axes):
    data = np.loadtxt(fn)
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    vix = np.round(data[:,0])
    vortex_idx = 0
    while np.any(vix == vortex_idx):
        ix = vix == vortex_idx

        plot_axes.plot(x[ix], y[ix], z[ix], '-o', ms=3)
        vortex_idx += 1
