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

import numpy as np
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
    parser.add_argument("--output-dir", help="Where to save picutes. data_dir by default")
    parser.add_argument("--slow", help="Plot segment-by-segment. Very slow, but right now needed for periodic boxes.",
                        action = "store_true");
    parser.add_argument("-D", help="Size of the plot box. The box will be interval [-D,D]^3. If the box is assymetric, use also -D1",
                       type = float)
    parser.add_argument("-D1", help="For assymetric plot boxes. The interval will be [D1, D]^3.",
                       type = float)
    parser.add_argument("--config", help="Path to the config file of the simulation. This will override D, D1 and slow.")
    parser.add_argument("--fix-plot-box", help="Use a plot dimensions fixed by the computation box in the config fig.",
                        action = "store_true")
    parser.add_argument('--show-time', help="Show time on the plots.",
                        action='store_true')
    parser.add_argument('--colorful', help='Plot vortices using many colors',
                        action='store_true')
    parser.add_argument('--projection', help='Projection axis.', type=str, default='z')
    parser.add_argument('--just-one', help='Plot just a single vortex.',
                        type=int, default=-1)
    parser.add_argument('--auto-ax', help='Leave axes range automatic',
                        action='store_true')
    parser.add_argument('--xlim', help='x-axis limiting values', type=float,
                        nargs=2)
    parser.add_argument('--ylim', help='y-axis limiting values', type=float,
                        nargs=2)
    parser.add_argument('--zlim', help='z-axis limiting values', type=float,
                        nargs=2)


    args = parser.parse_args()

    data_dir = args.data_dir
    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        output_dir = data_dir
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
        LBB = np.array(domain[0])
        RFT = np.array(domain[1])

        if args.fix_plot_box:
            Dxl, Dyl, Dzl = LBB*10
            Dxh, Dyh, Dzh = RFT*10
        else:
            Lmax = np.abs(LBB - RFT).max()
            mids = 0.5*(LBB + RFT)
            Dxl = mids[0] - Lmax/2
            Dxh = mids[0] + Lmax/2
            Dyl = mids[1] - Lmax/2
            Dyh = mids[1] + Lmax/2
            Dzl = mids[2] - Lmax/2
            Dzh = mids[2] + Lmax/2
        
        Dls = [Dxl, Dyl, Dzl]
        Dhs = [Dxh, Dyh, Dzh]
        dl_max = 2*config.dl_max
        
        dt = config.dt
        nshots = config.frame_shots
        
    #pick the plot axes based on the projection axis
    if args.projection == 'z':
        axids = (0, 1)
    elif args.projection == 'x':
        axids = (1, 2)
    elif args.projection == 'y':
        axids = (0,2)
    else:
        raise RuntimeError('Only x, y, z projections are possible.')


    files = glob(path.join(data_dir, 'frame*.dat'))
    files.sort()

    fig, ax = plt.subplots(1,1, figsize=(8,4.5), tight_layout=True)
        
    if args.colorful:
        color = None
    else:
        color = 'r'

    i=-1
    for fn in files:
        i = i+1
        dst_file = path.join(output_dir, path.split(fn)[-1].replace('.dat', '_2D.png'))
        if path.isfile(dst_file):
            continue
        print("{}/{}".format(i, len(files)))
        frame_index = int(path.split(fn)[1][5:-4])
        time = frame_index*nshots*dt

        ax.clear()
        if args.show_time:
            txt = fig.text(0.05, 0.05, "$t$ = {:.06f} s".format(time), fontsize=18)
        draw_vortices(fn, ax, slow=slow, max_len=dl_max, color=color, projection=axids)
        if args.xlim is None:
            ax.set_xlim(Dls[axids[0]], Dhs[axids[0]])
        else:
            ax.set_xlim(*args.xlim)
        if args.ylim is None:
            ax.set_ylim(Dls[axids[1]], Dhs[axids[1]])
        else:
            ax.set_ylim(*args.ylim)
        ax.set_aspect('equal')
        ax.set_xlabel("$x$ (mm)")
        ax.set_ylabel("$y$ (mm)")
        fig.savefig(dst_file)
        if args.show_time:
            txt.remove()
