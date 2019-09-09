#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>

This file is part of OpenVort.

OpenVort is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenVort is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenVort.  If not, see <http://www.gnu.org/licenses/>.

This file contains some basic utilities to work with the frame data.
"""

import numpy as np
import os.path as path

def frame_id(fn):
    fnb = path.split(fn)[-1]
    return int(fnb[5:-4])

def build_vortex(frame_data, vortex_idx, max_l = None):
    vidx = frame_data[:,0].astype(int)
    
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
        if next_node_i >= 0:
            k = np.where(node_i == next_node_i)[0][0]
            
        if next_node_i == initial_node_i:
            vxs.append(xs[k,:])
            break #we ran full circle
            
        if next_node_i < 0:
            k = 0
            while True:
                next_node_i = reverse[k]
                if next_node_i >= 0:
                    k = np.where(node_i == next_node_i)[0][0]
                    vxs = [xs[k]] + vxs
                else:
                    break
            break
    vxs = np.array(vxs)
    if max_l is not None:
        return split_vortex(vxs, max_l)
    return vxs

def split_vortex(vxs, max_l):
    xs = vxs[:,0]
    ys = vxs[:,1]
    zs = vxs[:,2]
    dx = np.diff(xs)
    dy = np.diff(ys)
    dz = np.diff(zs)
    
    d = np.sqrt(dx**2 + dy**2 + dz**2)
    jumps = np.nonzero(d > max_l)[0]
    if len(jumps) == 0:
        return [vxs]
    
    pieces = []
    k0 = 0
    for jump in jumps:
        pieces.append(vxs[k0:jump+1])
        k0 = jump+1
    pieces.append(vxs[k0:])
    
    return pieces

def stitch_vortex(pieces, box_Ds):
    """
    move pieces[1:] so that they form a continuous vortex starting with pieces[0]
    """
    full_vortex = pieces[0]
    for piece in pieces[1:]:
        vortex_delta = piece[0,:] - full_vortex[-1,:]
        
        #constuct the shift
        dir_flags = np.abs(vortex_delta) > box_Ds
        shift = -np.sign(vortex_delta)*box_Ds
        shift[np.logical_not(dir_flags)] = 0
        np.row_stack((full_vortex, piece + shift))
    
    return full_vortex