# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 16:37:22 2015

@author: emil
"""

import numpy as np
import scipy.interpolate as interp
import os

dn, fn = os.path.split(__file__)
data = np.loadtxt(os.path.join(dn, "rhosn_T.dat"))

T = data[:,0]
rhos = data[:,1]
rhon = data[:,2]

rs = interp.InterpolatedUnivariateSpline(T, rhos, k = 1)
rn = interp.InterpolatedUnivariateSpline(T, rhon, k = 1)
r  = lambda T: rs(T) + rn(T)