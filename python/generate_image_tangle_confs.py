#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 17:20:18 2018

@author: emil
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