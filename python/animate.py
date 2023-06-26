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
    parser.add_argument("--slow", help="Plot segment-by-segment. Very slow, but right now needed for periodic boxes.",
                        action = "store_true")
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
    parser.add_argument('--aspect', help='equal aspect ratio on the axes', action='store_true')

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
        LBB = np.array(domain[0])
        RFT = np.array(domain[1])

        if args.fix_plot_box:
            Dxl, Dyl, Dzl = LBB*10000
            Dxh, Dyh, Dzh = RFT*10000
        else:
            Lmax = np.abs(LBB - RFT).max()
            mids = 0.5*(LBB + RFT)
            Dxl = mids[0] - Lmax/2
            Dxh = mids[0] + Lmax/2
            Dyl = mids[1] - Lmax/2
            Dyh = mids[1] + Lmax/2
            Dzl = mids[2] - Lmax/2
            Dzh = mids[2] + Lmax/2
        dl_max = 2*config.dl_max
        
        dt = config.dt
        nshots = config.frame_shots


    files = glob(path.join(data_dir, 'frame*.dat'))
    files.sort()

    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='3d')
    
    if args.colorful:
        color = None
    else:
        color = 'r'
        
    if args.just_one >= 0:
        just_one = args.just_one
    else:
        just_one = None
        
    if args.xlim is not None:
        Dxl, Dxh = args.xlim
    if args.ylim is not None:
        Dyl, Dyh = args.ylim
    if args.zlim is not None:
        Dzl, Dzh = args.zlim

    i=-1
    for fn in files:
        i = i+1
        pic_path = fn.replace('.dat', '.png')
        if just_one is not None:
            pic_path = fn.replace('.dat', '_{:d}.png'.format(just_one))
        if path.isfile(pic_path):
            continue
        print("{}/{}".format(i, len(files)))
        frame_index = int(path.split(fn)[1][5:-4])
        time = frame_index*nshots*dt

        ax.clear()
        ax.auto_scale_xyz(1, 1, 1)
        draw_vortices(fn, ax, slow=slow, max_len=dl_max, color=color,
                      just_one=just_one)
        if not args.auto_ax:
            ax.set_xlim(Dxl, Dxh)
            ax.set_ylim(Dyl, Dyh)
            ax.set_zlim(Dzl, Dzh)
        ax.set_aspect('auto')
        ax.set_xlabel("$x$ (mm)")
        ax.set_ylabel("$y$ (mm)")
        ax.set_zlabel("$z$ (mm)")
        if args.aspect:
            ax.set_aspect('equal')
        if args.show_time:
            txt = fig.text(0.05, 0.05, "$t$ = {:.06f} s".format(time), fontsize=18)
        fig.tight_layout()
        fig.savefig(pic_path, dpi=120)
        if args.show_time:
            txt.remove()
