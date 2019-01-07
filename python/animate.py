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

import io, libconf

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Draws the individual frames.")
    parser.add_argument("data_dir", help="Directory with the frames.")
    parser.add_argument("--slow", help="Plot segment-by-segment. Very slow, but right now needed for periodic boxes.",
                        action = "store_true");
    parser.add_argument("-D", help="Size of the plot box. The box will be interval [-D,D]^3. If the box is assymetric, use also -D1",
                       type = float)
    parser.add_argument("-D1", help="For assymetric plot boxes. The interval will be [D1, D]^3.",
                       type = float)
    parser.add_argument("--config", help="Path to the config file of the simulation. This will override D, D1 and slow.")

    args = parser.parse_args()

    data_dir = args.data_dir
    slow = args.slow

    if args.D:
        D = args.D
    else:
        D = 1

    if args.D1:
        D1 = args.D1
    else:
        D1 = -D

    Dxl = D1; Dyl = D1; Dzl = D1
    Dxh = D; Dyh = D; Dzh = D
    dl_max = 0.05

    if args.config:
        with io.open(args.config) as f:
            config = libconf.load(f)
        domain = config['domain']
        LBB = domain[0]
        RFT = domain[1]

        Dxl = LBB[0]
        Dxh = RFT[0]
        Dyl = LBB[1]
        Dyh = RFT[1]
        Dzl = LBB[2]
        Dzh = RFT[2]
        dl_max = config.dl_max


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
        draw_vortices(fn, ax, slow=slow, max_len=dl_max)
        ax.set_xlim(Dxl, Dxh)
        ax.set_ylim(Dyl, Dyh)
        ax.set_zlim(Dzl, Dzh)
        ax.set_aspect('equal')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        fig.tight_layout()
        fig.savefig(fn.replace('.dat', '.png'), dpi=120)
