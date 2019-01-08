#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 11:18:01 2019

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
import argparse

def density_dist(data, z_max, res, D, dl_max):
    vidx = 0

    xs = np.linspace(-D, D, res)
    ys = np.linspace(-D, D, res)

    xs = 0.5*(xs[:-1] + xs[1:])
    ys = 0.5*(ys[:-1] + ys[1:])

    XX, YY = np.meshgrid(xs, ys)
    LS = np.zeros_like(XX)

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


        for (cx, cy, ls) in zip(cents[:,0], cents[:,1], lens):
            LS[np.argmin(np.abs(xs - cx)), np.argmin(np.abs(ys - cy))] += ls

        vidx += 1

    return XX, YY, LS

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Calculate the density profile of vortex line length along the z-direction.')
    parser.add_argument('filenames', help='The frame*.dat file with the tangle', nargs='+')
    parser.add_argument('--res', help='Resolution.', type=int, default = 10)
    parser.add_argument('--D', help='Domain half-size', type=float, default=0.05)
    parser.add_argument('--zmax', help='Highest layer.', type=float, default = 0.1)
    parser.add_argument('--dl_max', help='Maximum discretisation length.', type=float, default = 0.002)
    parser.add_argument('--plot', help='Plot the density profile.', action='store_true')
    parser.add_argument('--no-overwrite', help='Do not overwrite existing files.', action='store_true')

    args = parser.parse_args()

    XX_avg = 0
    YY_avg = 0
    LS_avg = 0

    for file in args.filenames:
        file_out = path.join(path.dirname(file), 'XY_dist'+path.split(file)[1])+'.npz'
        print(file_out)
        if path.isfile(file_out):
            print('Skipping ' + file_out)
            continue;
        data = np.loadtxt(file)
        XX, YY, LS = density_dist(data, args.zmax, args.res, args.D, args.dl_max)
        np.savez(file_out, XX=XX, YY=YY, LS=LS)
        XX_avg += XX
        YY_avg += YY
        LS_avg += LS

    XX_avg /= len(args.filenames)
    YY_avg /= len(args.filenames)
    LS_avg /= len(args.filenames)

    if args.plot:
        f, ax = plt.subplots(1,1)
        ax.contourf(XX_avg, YY_avg, LS_avg, 20)
#        ax.imshow(LS_avg)
        plt.show()