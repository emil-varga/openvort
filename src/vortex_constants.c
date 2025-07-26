/*******************************************************************************
 * Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>
 * 
 * This file is part of OpenVort
 * 
 * OpenVort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "vortex_constants.h"
#include "tangle.h"
//declared as extern in vortex_constants.h

int global_max_steps = -1;

double VORTEX_WIDTH = 1e-8; //cm
//quantum of circulation
double KAPPA = 9.97e-4; //cm^2/s

//Local induction approximation
int global_LIA_only = 0; //default false

//mutual friction parameters
double global_alpha = 0.1;
double global_alpha_p = 0.01;

//densities, g/cm^3
double global_rho_n = 1.0; //normal fluid
double global_rho_s = 0.45; //superfluid

//discretisation and stepping
double global_dt = 1e-3;
double global_dl_min = 1e-3;
double global_dl_max = 5e-3;
double global_reconnection_angle_cutoff = 0.087; //in radians, about 5 deg
double rec_dist = 1e-3;

int global_small_loop_cutoff = 5;
int global_frame_shot = 100;
int global_num_threads = 4;

//switch mutual friction on/off
int global_use_mutual_friction = 1;

int global_hyperfriction = 1;
double global_max_curvature_scale = 0.2;
double global_hyperalpha = 0.5;

/*
 * Pinning mode
 */
int global_pin_mode = PINNED;

//elimination of loops near the origin (for spherical flows)
int global_eliminate_origin_loops = 0; //default off
double global_eliminate_loops_origin_cutoff = 3e-2;

//similar to the spherical case, but for cylinder around zaxis
int global_eliminate_zaxis_loops = 0;
double global_eliminate_loops_zaxis_cutoff = 2e-2;

//remove loops further from the z-axis than the given cutoff
int global_eliminate_outer_loops = 0; //default off
double global_eliminate_outer_loops_cutoff = 1; //cm

/*
 * configuration of vortex injection
 */
//injecting loops at upper z-plane
int global_loop_injection = 0; //bool, inject or not
double global_loop_injection_frequency = 1; //injections per second

//line pair injection
int global_line_injection = 0; //bool, inject or not
int global_line_injection_n = 1; //how many pairs to inject
double global_line_injection_frequency; //injections per second
int global_line_injection_polarized = 0; //whether to inject the vortices in a polarized way, off by default

//Barnes-Hut approximation
int global_use_BH = 0; //false by default
int global_BH_quadtree = 0; //3D octree by default
double global_BH_resolution = 0.1; //maximum allowed ratio of BH box size to distance to center of mass
double global_BH_grain = 1e-8; //minimum size of the BH box
