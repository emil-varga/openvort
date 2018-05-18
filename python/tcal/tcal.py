# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 19:40:25 2015

@author: emil
"""

import numpy as np
import os

def get_cal_curve(n):
    dn, fn = os.path.split(__file__)
    d = np.loadtxt(os.path.join(dn, 'TR'+str(n)+'.cal'))
    n = int(d[0])
    p = d[1:n+2]
    s = d[-1]
    return (lambda R: 1/np.polyval(p, np.log(R))), s