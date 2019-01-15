#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 12:16:51 2019

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
import argparse

from density_profile_avg import get_frame_index

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Average the density xy-distribution calculated using density_xy_distribution.')
    parser.add_argument('filenames', help='Files to work on.', nargs='+')
    parser.add_argument('--start_ix', help='Starting index for averaging.', type=int, default=0)
    parser.add_argument('--end_ix', help='Last index for averaging. -1 for last', type=int, default=-1)
    parser.add_argument('--plot', help='Plot the result.', action='store_true')

    args = parser.parse_args()

    indices = [get_frame_index(file) for file in args.filenames]

    ix = np.argsort(indices)
    indices = np.array(indices)[ix]
    files = np.array(args.filenames)[ix]

    if args.end_ix >= 0:
        inrange_ix = np.logical_and(indices >= args.start_ix,
                            indices < args.end_ix)
    else:
        inrange_ix = indices >= args.start_ix

    XX = 0
    YY = 0
    LS = 0
    LS2 = 0
    m = len(files[inrange_ix])
    for file in files[inrange_ix]:
        data = np.load(file)
        XX += data['XX']
        YY += data['YY']
        LS += data['LS']
        LS2 += data['LS']**2

    XX /= m
    YY /= m
    LS /= m
    LS2 = np.sqrt(LS2/m - LS**2)

    dirname = path.dirname(args.filenames[0])
    file_out = 'Lxydist_avg_{}_{}.npz'.format(args.start_ix, args.end_ix)
    np.savez(path.join(dirname, file_out), XX = XX, YY = YY, LS=LS, LS2=LS2)

    if args.plot:
        f, ax = plt.subplots(1,1)
#        ax.contourf(XX, YY, LS, 200)
        c = ax.imshow(LS, interpolation='bilinear', origin='lower',
                  extent=(XX.min(), XX.max(), YY.min(), YY.max()))
        plt.colorbar(c)
        ax.set_aspect('equal')
        plt.show()