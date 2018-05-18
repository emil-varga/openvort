#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 09:35:40 2017

@author: emil
"""

import numpy as np
import scipy.interpolate as intp

import os

dn, fn = os.path.split(__file__)
data = np.loadtxt(os.path.join(dn, "schwarz_coefs.dat"))

T = data[:,0]
_c1 = data[:,1]
_c2 = data[:,2]
_Il = data[:,3]
_alpha = data[:,4]

c1 = intp.InterpolatedUnivariateSpline(T, _c1, k = 3)
c2 = intp.InterpolatedUnivariateSpline(T, _c2, k = 3)
Il = intp.InterpolatedUnivariateSpline(T, _Il, k = 3)

def alpha(Ts):
    xalpha = intp.InterpolatedUnivariateSpline(T, np.log(_alpha), k = 3)
    return np.exp(xalpha(Ts))

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.close('all')
    Ts = np.linspace(1, 2.15, 1000)
    alphas = alpha(Ts)
    
    f, ax = plt.subplots(1,1)
    ax.semilogx(alphas, Ts)
    ax.semilogx(_alpha, T, 'o')
    
    c2s = c2(Ts)
    f, ax = plt.subplots(1,1)
    ax.plot(Ts, c2s)
    ax.plot(T, _c2, 'o')