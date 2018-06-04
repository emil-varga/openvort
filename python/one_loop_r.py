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

import os.path as path
from glob import glob

plt.close('all')

data_dir = '../data_oneloop'
files = glob(path.join(data_dir, 'frame*.dat'))
files.sort()

vs = []
rs = []
rse = []
rcs = []
ts = []
vs_teor = []
alphap = 0.02
#alphap=0

a = 1e-8
kappa = 9.97e-4;

for fn in files:
    d = np.loadtxt(fn)

    x = d[:,1]
    y = d[:,2]
    z = d[:,3]

    vx = d[:,4]
    vy = d[:,5]
    vz = d[:,6]

    x0 = x.mean()
    y0 = y.mean()
    z0 = z.mean()

    ts.append(abs(d[:,7]*d[:,10] + d[:,8]*d[:,11] + d[:,9]*d[:,12]).max())

    r_mean = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2).mean()
    r_dev =  np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2).std()
    r_crv_mean = 1/np.sqrt(d[:,10]**2 + d[:,11]**2 + d[:,12]**2).mean()
#    vx_mean = np.sqrt(vx**2 + vy**2 + vz**2).mean()
    vx_mean = vx.mean()

    v_teor = kappa/4/np.pi/r_mean*(np.log(8*r_mean/a) - 0.5)*(1-alphap)

    rs.append(r_mean)
    vs.append(vx_mean)
    vs_teor.append(v_teor)
    rcs.append(r_crv_mean)
    rse.append(r_dev)


f, ax = plt.subplots(1,1)
ax.plot(rs, vs, 'o', label = 'numerics, alpha = 0.1, alphap = 0.2')
ax.plot(rs, vs_teor, '-', lw=3, label='mf-corrected theory, alphap = 0.2')
ax.plot(rs, np.array(vs_teor)/(1-alphap), '--', lw=3, label='theory, no mf')
ax.legend(loc='best')
ax.set_xlabel('radius (cm)')
ax.set_ylabel('velocity (cm/s)')
f.tight_layout()

f, ax = plt.subplots(1,1)
vs = np.array(vs)
vs_teor = np.array(vs_teor)
ax.plot(rs, 2*100*abs(vs - vs_teor)/(vs + vs_teor), '.')