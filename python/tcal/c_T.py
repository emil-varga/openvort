#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 13:11:31 2017

@author: emil
"""

import numpy as np
import scipy.interpolate as interp
import os

dn, fn = os.path.split(__file__)
data = np.loadtxt(os.path.join(dn, "T_c.dat"))

T = data[:,0]
c = data[:,1] * 0.249837 #units J/g/K

c = interp.InterpolatedUnivariateSpline(T, c, k = 1)