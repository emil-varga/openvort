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

def draw_vortices(fn, plot_axes, slow=False, max_len=0.05):
    data = np.loadtxt(fn)
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    vix = np.round(data[:,0])
    vortex_idx = 0
    while np.any(vix == vortex_idx):
        ix = vix == vortex_idx
        if slow:
            rs = np.column_stack((x[ix], y[ix], z[ix]))
            xs = []
            ys = []
            zs = []
            for r, rn in zip(rs[:-1,:], rs[1:,:]):
                if np.linalg.norm(r - rn) < max_len:
                    xs.append(r[0])
                    ys.append(r[1])
                    zs.append(r[2])
                else:
                    plot_axes.plot(xs, ys, zs, '-', color='red', lw=0.5)
                    xs, ys, zs = [], [], []
            plot_axes.plot(xs, ys, zs, '-', color='red', lw=0.5)
            if np.linalg.norm(rs[0,:] - rs[-1,:]) < max_len:
                plot_axes.plot([rs[0,0], rs[-1,0]],
                               [rs[0,1], rs[-1,1]],
                               [rs[0,2], rs[-1,2]],
                               '-', color='red', lw = 0.5)
        else:
            plot_axes.plot(x[ix], y[ix], z[ix], '-', ms=2, lw=0.5,
                           color = 'red')
        vortex_idx += 1

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import argparse
    parser = argparse.ArgumentParser('Utility to plot the vortex tangle.')
    parser.add_argument('filename', help='Filename of the .dat file with the vortices.')
    parser.add_argument('--dl_max', help='Maximum length between vortex points.', default=0.05,
                        type=float)

    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    draw_vortices(args.filename, ax, slow=True, max_len=args.dl_max)
    plt.show()