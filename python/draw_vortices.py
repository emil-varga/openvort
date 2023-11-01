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

def draw_vortices(fn, plot_axes, slow=False, max_len=0.05, scale=10, just_one=None,
                  color=None, colorcode_z=True, projection=None, **pltkw):
    data = np.loadtxt(fn)

    vix = data[:,0].astype(int)
    vortex_idx = 0
    if projection is None:
        plot_axes.set_proj_type('ortho')
    if just_one is not None:
        vortex_idx = just_one
    while np.any(vix == vortex_idx):
        vxs_pieces = util.build_vortex(data, vortex_idx, max_l=max_len)
        clr = color
        for vxs in vxs_pieces:
            vxs *= scale
            dd = vxs.shape
            if dd[0] < 2:
                continue
            if colorcode_z:
                if vxs[dd[0]-1,2] > vxs[0,2]:
                    clr = 'r'
                elif vxs[dd[0]-1,2] < vxs[0,2]:
                    clr = 'b'
                else:
                    clr = 'g' 
            if projection is None:
                pl = plot_axes.plot(1000*vxs[:,0], 1000*vxs[:,1], 1000*vxs[:,2], '-', color=clr, **pltkw)
            else:
                ax1, ax2 = projection
                pl = plot_axes.plot(vxs[:,ax1], vxs[:,ax2], '-', color=clr, lw=2)
                plot_axes.plot([vxs[:,ax1].mean()], [vxs[:,ax2].mean()], 'x', color=clr, ms=2, **pltkw)
            if clr is None:
                clr = pl[-1].get_color()
        if just_one is not None:
            break
        vortex_idx += 1

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import argparse
    parser = argparse.ArgumentParser('Utility to plot the vortex tangle.')
    parser.add_argument('filename', help='Filename of the .dat file with the vortices.')
    parser.add_argument('--dl_max', help='Maximum length between vortex points.', default=0.05,
                        type=float)
    parser.add_argument('--just-one', help='Plot just a single vortex.',
                        type=int, default=-1)

    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if args.just_one >= 0:
        just_one = args.just_one
    else:
        just_one = None
    draw_vortices(args.filename, ax, slow=False, max_len=args.dl_max, just_one=just_one)
    ax.set_xlabel('$x$ ($\mu$m)')
    ax.set_ylabel('$y$ ($\mu$m)')
    ax.set_zlabel('$z$ ($\mu$m)')
    
    fig.tight_layout()
    fig.savefig(args.filename.replace('.dat', '.pdf'))
    plt.show()
