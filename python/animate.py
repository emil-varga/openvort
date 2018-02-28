#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 13:50:04 2017

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os.path as path
from glob import glob

from draw_vortices import draw_vortices

data_dir = '../data_rec'

files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

i=-1
for fn in files:
    i = i+1
    if path.isfile(fn.replace('.dat', '.png')):
        continue
    print("{}/{}".format(i, len(files)))

    ax.clear()
    draw_vortices(fn, ax)
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(-0.25, 0.25)
    ax.set_zlim(-0.25, 0.25)
    ax.set_aspect('equal')
    fig.savefig(fn.replace('.dat', '.png'))

#plt.close(fig)
