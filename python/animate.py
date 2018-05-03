#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 13:50:04 2017

@author: emil
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os.path as path
from glob import glob

from draw_vortices import draw_vortices
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('provide data dir as the command line argument')
    data_dir = path.abspath(sys.argv[1])
    
    slow=False
    if len(sys.argv) > 2:
        if sys.argv[2] == "slow":
            slow = True

    files = glob(path.join(data_dir, 'frame*.dat'))
    files.sort()

    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='3d')

    i=-1
    for fn in files:
        i = i+1
        if path.isfile(fn.replace('.dat', '.png')):
            continue
        print("{}/{}".format(i, len(files)))

        ax.clear()
        ax.auto_scale_xyz(1, 1, 1)
        draw_vortices(fn, ax, slow=slow)
        ax.set_xlim(-0.25, 0.25)
        ax.set_ylim(-0.25, 0.25)
        ax.set_zlim(-0.25, 0.25)
        ax.set_aspect('equal')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        fig.tight_layout()
        fig.savefig(fn.replace('.dat', '.png'), dpi=120)
