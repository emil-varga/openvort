#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 18:07:11 2018

@author: emil
"""

config_template = """
domain = ([-0.05, -0.05, -0.05], [0.05, 0.05, 0.05]);

boundaries = "open";

use_mutual_friction = True;
alpha = 0.034;
alpha_p = 1.383e-2;

frame_shots = 1;
dt = 1e-4; #seconds
dl_min = 1e-3; #cm
dl_max = 2e-3; #cm
small_loop_cutoff = 5; #points
reconnection_angle_cutoff = 5; #degrees
reconnection_distance = 2e-3;
num_threads = 4
eliminate_origin_loops = True;
eliminate_loops_origin_cutoff = 2e-2

init_mode = "restart";
init_file = "T130/data_1mms@5mm_olr_fix4/frame6000.dat"

#vns = 0.1 cm/s 5 mm away from the origin

vn_conf = {
        type     = "no flow";
        strength = 0.300044613;
        cutoff   = 2e-2;
};

vs_conf = {
        type     = "no flow";
        strength = -0.014114653;
        cutoff   = 2e-2;
};
"""

