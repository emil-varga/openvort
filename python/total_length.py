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

import sys
import libconf

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('provide data dir as the command line argument')
    data_dir = path.abspath(sys.argv[1])
    cfg_filename = glob(path.join(data_dir, '*.cfg'))[0]
    
    with open(path.join(data_dir, cfg_filename)) as cfg_file:
        config = libconf.load(cfg_file);
        
    dt0 = config['dt']
    frame_shots = config['frame_shots']
    dt = frame_shots * dt0
    
    files = glob(path.join(data_dir, 'frame*.dat'))
    files.sort()
    
    time = dt
    times = []
    lengths = []
    
    for fn in files:
        print(fn)
        d = np.loadtxt(fn)
        vidx = 0
        l = 0
        while True:
            ix = d[:,0] == vidx;
            if not any(ix):
                break
    
            rs = d[ix, 1:4]
            dr = np.diff(rs, axis=0)
            l += np.sum(np.sqrt(dr**2))
            vidx += 1
        lengths.append(l)
        times.append(time)
        time += dt
    
    out = np.column_stack((times, lengths))
    np.savetxt(path.join(data_dir, 'L_t.txt'), out, header = 'time(s)\ttotal length(cm)')