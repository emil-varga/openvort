# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:14:11 2015

@author: emil
"""

from .tcal import get_cal_curve
from .tofp import get_T_interpolant
from .rhosnoft import rs, rn, r
from .soft import s
from .c_T import c
from .mutual_friction import B, Bp, alpha, alphap

from .schwarz_coefs import c1, c2, Il, alpha

tp = get_T_interpolant()

def vn(A, T, Q):
    return Q/(A*s(T)*T*r(T))

def vs(A, T, Q):
    return vn(A, T, Q)*rn(T)/rs(T)

def vns(A, T, Q):
    return abs(vn(A, T, Q)) + abs(vs(A, T, Q))