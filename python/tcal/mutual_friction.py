# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:15:41 2016

@author: emil
"""

import numpy as np
import scipy.interpolate as interp
import os

dn, fn = os.path.split(__file__)
data = np.loadtxt(os.path.join(dn, "mutual_friction.dat"))

T_data = data[:,0]
B_data = data[:,1]
Bp_data = data[:,2]
alpha_data = data[:,3]
alphap_data = data[:,4]

B = interp.InterpolatedUnivariateSpline(T_data, B_data, k = 1)
Bp = interp.InterpolatedUnivariateSpline(T_data, B_data, k = 1)
alpha = interp.InterpolatedUnivariateSpline(T_data, B_data, k = 1)
alphap = interp.InterpolatedUnivariateSpline(T_data, B_data, k = 1)