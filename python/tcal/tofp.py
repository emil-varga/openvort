# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 19:01:03 2015

@author: emil
"""

import numpy as np
import scipy.interpolate as interp
import os

def torr_to_Pa(torrP):
    return torrP*101325.0/760.0

def get_T_interpolant():
    dn, fn = os.path.split(__file__)
    data = np.loadtxt(os.path.join(dn, "T_P.dat"))
    
    T = data[:,0]
    P = data[:,1]
    
    Tp = interp.InterpolatedUnivariateSpline(P, T, k = 1)
    
    return (lambda P: Tp(torr_to_Pa(P)))