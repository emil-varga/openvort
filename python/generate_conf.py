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
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

The purpose of these functions is to generate config file for the OpenVort run.
This is useful especially for the mutual friction parameters which can be looked up
automatically using the temperature from tcal module.
"""

import tcal

config_template_general = """
domain = ([-{D}, -{D}, -{D}], [{D}, {D}, {D}]);

boundaries = "{boundary}";

use_mutual_friction = {use_mf};
alpha = {alpha};
alpha_p = {alphap};

frame_shots = {frame_shots};
dt = {dt}; #seconds
dl_min = {dl_min}; #cm
dl_max = {dl_max}; #cm
small_loop_cutoff = {small_cutoff}; #points
reconnection_angle_cutoff = {rec_angle}; #degrees
reconnection_distance = {rec_d};
num_threads = {nthreads}
eliminate_origin_loops = {origin_removal};
eliminate_loops_origin_cutoff = {origin_removal_cutoff};
"""

def generate_v_conf(vnvs, flow_type, parameters):
    v_conf = "{}_conf = {{\n\ttype = \"{}\";".format(vnvs, flow_type)

    for p in parameters:
        v_conf += '\n\t{} = {};'.format(p, parameters[p])
    v_conf += '\n};'

    return v_conf

def generate_init_conf(init_type, init_parameters):
    init_conf = "init_mode = {};".format(init_type)

    for p in init_parameters:
        init_conf += "\n{} = {};".format(p, init_parameters[p])

    return init_conf

def generate_general_conf(D = 0.05, boundary = 'open', use_mf = True, alpha = 0.1, alphap = 0.01,
                          frame_shots = 5, dt = 1e-4, dl_min = 1e-3, dl_max = 2e-3,
                          small_cutoff = 5, rec_angle = 5, rec_d = 2e-3,
                          nthreads = 4,
                          eliminate_origin_loops = False, origin_removal_cutoff = 2e-2):

    return config_template_general.format(D=0.1)

print(generate_v_conf('vn', 'spherical', {'str':10, 'cutoff':5}))