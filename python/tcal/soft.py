# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 16:43:32 2015

@author: emil
"""

import numpy as np
import scipy.interpolate as interp
import os

dn, fn = os.path.split(__file__)
data = np.loadtxt(os.path.join(dn, "s_T.dat"))

T = data[:,0]
ss = data[:,1]

s = interp.InterpolatedUnivariateSpline(T, ss, k = 1)