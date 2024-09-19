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
from glob import glob
import os.path as path

from util import frame_id

import argparse

import libconf

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Calculates the total lengt of the vortices.')
    parser.add_argument('data_dir', help='Directory with the frames.')
    parser.add_argument('--config', help='Config file.')
    args = parser.parse_args()
    data_dir = args.data_dir
    cfg_filename = args.config

    with open(cfg_filename) as cfg_file:
        config = libconf.load(cfg_file);

    dt0 = config['dt']
    frame_shots = config['frame_shots']
    dt = frame_shots * dt0

    files = glob(path.join(data_dir, 'frame*.dat'))
    files.sort(key=frame_id)

    time = dt
    times = []
    lengths = []

    starti = 0
    output_file = path.join(data_dir, 'L_t.txt')
    if path.isfile(output_file):
        Lt = np.loadtxt(output_file)
        starti = Lt.shape[0] #1 line per file
        time = (starti+1)*dt

    for fn in files[starti:]:
        print(fn)
        d = np.loadtxt(fn)
        vidx = 0
        fwd = d[:,-4]
        l = np.sum(d[fwd>=0,-1])
        lengths.append(l)
        times.append(time)
        time += dt

    out = np.column_stack((times, lengths))
    
    if path.isfile(output_file):
        out = np.row_stack((Lt, out))

    f, ax = plt.subplots(1,1)
    ax.plot(out[:,0], out[:,1])
    ax.set_xlabel('time (s)')
    ax.set_ylabel('total length (cm)')
    plt.show()

    np.savetxt(path.join(data_dir, 'L_t.txt'), out, header = 'time(s)\ttotal length(cm)')
