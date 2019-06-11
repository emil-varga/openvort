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
import util

def draw_vortices(fn, plot_axes, slow=False, max_len=0.05, scale=10, just_one=False):
    data = np.loadtxt(fn)

    vix = np.round(data[:,0])
    vortex_idx = 0
    plot_axes.set_proj_type('ortho')
    if just_one: vortex_idx = np.random.choice(int(vix.max()))
    while np.any(vix == vortex_idx):
        vxs_pieces = util.build_vortex(data, vortex_idx, max_l=max_len)
        for vxs in vxs_pieces:
            vxs *= scale
            plot_axes.plot(vxs[:,0], vxs[:,1], vxs[:,2], '-', color='r')        
        if just_one: break
        vortex_idx += 1

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import argparse
    parser = argparse.ArgumentParser('Utility to plot the vortex tangle.')
    parser.add_argument('filename', help='Filename of the .dat file with the vortices.')
    parser.add_argument('--dl_max', help='Maximum length between vortex points.', default=0.05,
                        type=float)
    parser.add_argument('--just-one', help='Plot just a single random vortex.',
                        action = 'store_true')

    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    draw_vortices(args.filename, ax, slow=False, max_len=args.dl_max, just_one=args.just_one)
    plt.show()