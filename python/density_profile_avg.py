#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 10:31:15 2019

@author: emil
"""

import numpy as np
import matplotlib.pyplot as plt

import os.path as path
import argparse
import re

def get_frame_index(filename):
    fn = path.split(filename)[1]

    m = re.search('frame[0-9]{4}.dat', fn)
    if m is None:
        raise RuntimeError('Cannot decode filename ' + fn)

    group = m.group()
    index = int(group[5:9])

    return index

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Average the density profiles calculated using density_profile.')
    parser.add_argument('filenames', help='Files to work on.', nargs='+')
    parser.add_argument('--start_ix', help='Starting index for averaging.', type=int, default=0)
    parser.add_argument('--end_ix', help='Last index for averaging. -1 for last', type=int, default=-1)
    parser.add_argument('--plot', help='Plot the result.', action='store_true')
    parser.add_argument('--output', '-o', help='Output file basename.', default='Lprof')

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

    ZS = 0
    LS = 0
    LS2 = 0
    m = len(files[inrange_ix])
    for file in files[inrange_ix]:
        data = np.loadtxt(file)
        ZS += data[:,0]
        LS += data[:,1]
        LS2 += data[:,1]**2

    ZS /= m
    LS /= m
    LS2 = np.sqrt(LS2/m - LS**2)

    dirname = path.dirname(args.filenames[0])
    file_out = '{}_avg_{}_{}.dat'.format(args.output, args.start_ix, args.end_ix)
    np.savetxt(path.join(dirname, file_out),
               np.column_stack((ZS, LS, LS2)))

    if args.plot:
        f, ax = plt.subplots(1,1)
        ax.fill_between(ZS, LS-LS2, LS+LS2, alpha=0.5)
        ax.plot(ZS, LS, '-o')
        plt.show()