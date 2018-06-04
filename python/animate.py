#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>

This file is part of OpenVort.

OpenVort is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenVort is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenVort.  If not, see <http://www.gnu.org/licenses/>.
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
