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

#periodic
total = 0
for k in [-1, 0, 1]:
    for l in [-1, 0, 1]:
        for m in [-1, 0, 1]:
            shifts = np.sum(np.abs([k, l, m]))
            if shifts > 0 and shifts < 4: # <2 for 6, <3 for 18 <4 for 26 etc..
                print("{{{{{}, {}, {}}}, -1}},".format(k, l, m))
                total += 1
print(total)

#1 wall in Z-dir
total = 0
for k in [-1, 0, 1]:
    for l in [-1, 0, 1]:
        for m in [-1, 0, 1]:
            shifts = np.sum(np.abs([k, l, m]))
            if shifts > 0 and shifts < 4: # <2 for 6, <3 for 18 <4 for 26 etc..
                if m < 0:
                    mirror = 'Z_L'
                else:
                    mirror = '-1'
                print("{{{{{}, {}, {}}}, {}}},".format(k, l, m, mirror))
                total += 1
print(total)