#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 13:50:04 2017

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
from glob import glob

from draw_vortices import draw_vortices

data_dir = '../data'

files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for fn in files:
    if path.isfile(fn.replace('.dat', '.png')):
        continue
    
    ax.clear()
    draw_vortices(fn, ax)
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(-0.1, 0.1)
    ax.set_zlim(-0.5, 0.1)
    ax.set_aspect('equal')
    fig.savefig(fn.replace('.dat', '.png'))

plt.close(fig)