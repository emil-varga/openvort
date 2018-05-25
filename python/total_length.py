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
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os.path as path

data_dir = '../data_spherical'

files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

lengths = []
files = files[:]

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

f, ax = plt.subplots(1,1)
ax.plot(lengths)