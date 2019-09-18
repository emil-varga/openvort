#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (C) 2019 Emil Varga <varga.emil@gmail.com>

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

This script calculates the 2D vortex correlation index of a 2-dimensionalised tangle. The
tangle is reduced to 2 dimensions by averaging out the 'z' dimension.
"""

import numpy as np
import matplotlib.pyplot as plt
import util
import scipy.spatial as spat

import os.path as path
from glob import glob

import argparse
import libconf

def build_2d_projections(vortex_data, box_Ds = [1e-3, 1e-3, 1e-4]):
    """
    
    """
    node_ids = vortex_data[:,0].astype(int)
    
    vortex_ids = range(int(node_ids.min()), int(node_ids.max()))
    
    xs     = [] #vortex centre of mass, x component
    xss    = [] #vortex spread in x, standard deviation 
    ys     = [] #same as for y
    yss    = []
    kappas = [] #direction of circulation, +/- 1
    for vid in vortex_ids:
        vort_pieces = util.build_vortex(vortex_data, vid, max_l = 1e-5)
        vortex = util.stitch_vortex(vort_pieces, box_Ds = box_Ds) #10um x 10um x 1 um
        
        if vortex.shape[0] < 2:
            continue
        
        dz = vortex[-1,2] - vortex[0,2]
        if dz > box_Ds[2]*3/4:
            kappa = +1
        elif -dz > box_Ds[2]*3/4:
            kappa = -1
        else:
            #this is a vortex loop pinned at one wall
            continue
        
        kappas.append(kappa)
        xs.append(vortex[:,0].mean())
        xss.append(vortex[:,0].std())
        ys.append(vortex[:,1].mean())
        yss.append(vortex[:,1].std())
    
    xs     = np.array(xs)
    xss    = np.array(xss)
    ys     = np.array(ys)
    yss    = np.array(yss)
    kappas = np.array(kappas)
    
    return kappas, xs, ys, xss, yss

def calculate_2D_correlation_index(vortex_data, box_Ds, full_output=False):
    kappas, xs, ys, xss, yss = build_2d_projections(vortex_data, box_Ds)
    #KDTree expects the periodic domain to be [0, L]
    #so we shift the points to this box
    Lx = box_Ds[0]
    Ly = box_Ds[1]
    points = np.column_stack((xs + Lx/2, ys + Ly/2)) 
    #find the position of all 1st nearest neighbours
    tree = spat.cKDTree(points, boxsize = (Lx, Ly))
    rs, ids = tree.query(points, k=2)
    ids = ids[:,1]
#    print(rs.shape)
    Ci = np.sum(kappas*kappas[ids])/len(kappas)
    if full_output:
        return Ci, ids, kappas, xs, ys
    return Ci

if __name__ == '__main__':
    plt.close('all')
    
    parser = argparse.ArgumentParser(
            """Calculates the 2D correlation index for the frame*.dat files in 
            a given directory.""")
    parser.add_argument('directory', help='Directory with the frame*.dat data')
    parser.add_argument('--config', help='Config file.', required=True)
    parser.add_argument('--every', help='Work only on files[::every]. For speed.',
                        type=int, default=50)
    parser.add_argument('--append', help='Do not recalculate everything. Append to existing results file.',
                        action = 'store_true')
    parser.add_argument('--plot-nn', help='Plot the nearest neighbours for every frame.',
                        action = 'store_true')
    parser.add_argument('--no-save', help='Do not save the resulting Ci(t).',
                        action = 'store_true')
    
    args = parser.parse_args()
    directory   = args.directory
    config_file = args.config
    every       = args.every
    
    #pull the necessary info out of the config file
    print("Loading conf from {}".format(config_file))
    with open(config_file, 'r') as conf:
        config = libconf.load(conf)
    dt_step = config['dt'] #dt between simulation steps
    frame_shots = config['frame_shots']
    dt = dt_step * frame_shots #dt between succsessive frames
    
    domain = config['domain']
    Lx = domain[1][0] - domain[0][0]
    Ly = domain[1][1] - domain[0][1]
    Lz = domain[1][2] - domain[0][2]
    
    #find and sort the frame files
    files = glob(path.join(directory, 'frame*.dat'))
    files.sort(key=util.frame_id)
    
    #load the previous results, if needed
    if args.append:
        ts_all = np.array([dt*util.frame_id(file) for file in files])
        old_Cis = np.loadtxt(path.join(directory, 'Ci_t.txt'))
        ts = list(old_Cis[:,0])
        Cis = list(old_Cis[:,1])
        idx0 = np.where(ts_all > ts[-1])[0][0]
        files = files[idx0+1:]
    else:
        Cis = []
        ts = []
    
    files = files[::every]
    #start calculating the correlation index
    for k, file in enumerate(files):
        print("Progress: {}/{} ({})".format(k, len(files), path.split(file)[1]))
        data = np.loadtxt(file)
        frame_id = util.frame_id(file)
        ts.append(frame_id*dt)        
        if args.plot_nn:
            Ci, ids, kappas, xs, ys = calculate_2D_correlation_index(data, (Lx,Ly,Lz),
                                                                     full_output=True)
        else:
            Ci = calculate_2D_correlation_index(data, (Lx,Ly,Lz))
        Cis.append(Ci)
        
        if args.plot_nn:
            fig, ax = plt.subplots(1,1)
            ax.plot(xs[kappas>0], ys[kappas>0], 'ro')
            ax.plot(xs[kappas<0], ys[kappas<0], 'bo')
            ax.set_aspect('equal')
            
            for k, nk in enumerate(ids):
                x0 = xs[k]
                y0 = ys[k]
                x1 = xs[nk]
                y1 = ys[nk]
                ax.annotate("", xy=(x1, y1), xytext=(x0,y0), arrowprops=dict(arrowstyle="->"))
    fig, ax = plt.subplots(1,1)
    ax.plot(ts, Cis, '-')
    ax.axhline(0, color='k')
    if not args.no_save:
        out = np.column_stack((ts, Cis))
        outfile = path.join(directory, "Ci_t.txt")
        print(outfile)
        np.savetxt(outfile, out, header='time(s)\tCi')
    plt.show()