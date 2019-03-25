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

import argparse
import os.path as path

from glob import glob

def density_profile(data, layer_height, z_min, z_max, dl_max, window=None):
    vidx = 0

    ZS = np.arange(z_min, z_max, layer_height)
    LS = np.zeros_like(ZS)

    while True:
        ix = (data[:,0] == vidx)
        if not any(ix):
            break

        rvs = data[ix,1:4]

        cents = 0.5*(rvs[1:,:] + rvs[:-1,:])
        diffs = rvs[1:,:] - rvs[:-1, :]
        lens = np.sqrt(np.sum(diffs**2, axis=1))
        cents_z = cents[:,2]


        valid_ix = np.logical_and(lens < dl_max, cents_z < z_max);
        cents = cents[valid_ix,:]
        lens = lens[valid_ix]

        if window is not None:
            wx, wy, Dx, Dy = window
            cents_x = cents[:,0]
            cents_y = cents[:,1]
            valid_ix = np.logical_and(np.abs(cents_x - wx) < Dx/2,
                                      np.abs(cents_y - wy) < Dy/2)
            cents = cents[valid_ix,:]
            lens = lens[valid_ix]


        for (cz, ls) in zip(cents_z, lens):
            LS[np.argmin(np.abs(ZS - cz))] += ls

        vidx += 1

    return ZS, LS

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Calculate the density profile of vortex line length along the z-direction.')
    parser.add_argument('filenames', help='The frame*.dat file with the tangle', nargs='+')
    parser.add_argument('--dL', help='Layer height.', type=float, default = 0.01)
    parser.add_argument('--zmin', help='Lowest layer.', type=float, default = 0)
    parser.add_argument('--zmax', help='Highest layer.', type=float, default = 0.1)
    parser.add_argument('--dl_max', help='Maximum discretisation length.', type=float, default = 0.002)
    parser.add_argument('--plot', help='Plot the density profile.', action='store_true')
    parser.add_argument('--no-overwrite', help='Do not overwrite existing files.', action='store_true')
    parser.add_argument('--window', nargs=4, type=float,
                        help="Calculate the profile only in a square window in the xy coordinates.")

    args = parser.parse_args()

    ZS_avg = 0
    LS_avg = 0

    for file in args.filenames:
        dirname, basename = path.split(file)
        if args.window is not None:
            output_fn = path.join(dirname, 'WLprof_'+str(args.window)+'_'+basename)
        else:
            output_fn = path.join(dirname, 'Lprof_'+basename)
        if args.no_overwrite and path.exists(output_fn):
            continue
        print(output_fn)
        data = np.loadtxt(file)
        ZS, LS = density_profile(data, args.dL, args.zmin, args.zmax, args.dl_max, window=args.window)
        out = np.column_stack((ZS, LS))
        np.savetxt(output_fn, out)
        ZS_avg += ZS
        LS_avg += LS

    ZS_avg /= len(args.filenames)
    LS_avg /= len(args.filenames)

    if args.plot:
        f, ax = plt.subplots(1,1)
        ax.plot(ZS_avg, LS_avg, '-o')
        plt.show()