#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 15:57:22 2019

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os.path as path

from glob import glob

def density_profile(data, layer_height, z_min, z_max, dl_max):
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

    args = parser.parse_args()

    ZS_avg = 0
    LS_avg = 0

    for file in args.filenames:
        print(file)
        data = np.loadtxt(file)
        ZS, LS = density_profile(data, args.dL, args.zmin, args.zmax, args.dl_max)
        out = np.column_stack((ZS, LS))
        np.savetxt('Lprof_'+file, out)
        ZS_avg += ZS
        LS_avg += LS

    ZS_avg /= len(args.filenames)
    LS_avg /= len(args.filenames)

    if args.plot:
        f, ax = plt.subplots(1,1)
        ax.plot(ZS_avg, LS_avg, '-o')
        plt.show()