#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:29:48 2019

@author: emil
"""

import numpy as np

def build_vortex(frame_data, vortex_idx):
    vidx = frame_data[:,0]
    
    vortex_data = frame_data[vidx==vortex_idx,:]
    
    xs = vortex_data[:,1:4]
    
#    pin_wall = vortex_data[:,-1]
#    status_flag = vortex_data[:,-2]
    forward = vortex_data[:,-3]
    reverse = vortex_data[:,-4]
    node_i = vortex_data[:,-5]
    
    vxs = []
    
    k = 0
    initial_node_i = node_i[k]
    while True:
        vxs.append(xs[k,:])
        next_node_i = forward[k]
        if next_node_i == initial_node_i:
            break #we ran full circle
        if next_node_i >= 0:
            k = np.where(node_i == next_node_i)[0][0]
        else:
            k = 0
            while True:
                next_node_i = reverse[k]
                if next_node_i >= 0:
                    k = np.where(node_i == next_node_i)[0][0]
                    vxs = [xs[k]] + vxs
                else:
                    break
            break
    
    return np.array(vxs)        