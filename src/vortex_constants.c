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

double VORTEX_WIDTH = 1e-8; //cm
//quantum of circulation
double KAPPA = 9.97e-4; //cm^2/s

//mutual friction parameters
double alpha = 0.1;
double alpha_p = 0.01;

//densities, g/cm^3
double rho_n = 1.0; //normal fluid
double rho_s = 0.45; //superfluid

//discretisation and stepping
double global_dt = 1e-3;
double global_dl_min = 1e-3;
double global_dl_max = 5e-3;
double reconnection_angle_cutoff = 0.087; //in radians, about 5 deg
double rec_dist = 1e-3;

int small_loop_cutoff = 5;
int frame_shot = 100;
int global_num_threads = 4;

//switch mutual friction on/off
int use_mutual_friction = 1;

/*
 * Pinning mode
 */
int pin_mode = PINNED;

//elimination of loops near the origin (for spherical flows)
int eliminate_origin_loops = 0; //default off
double eliminate_loops_origin_cutoff = 3e-2;

//similar to the spherical case, but for cylinder around zaxis
int eliminate_zaxis_loops = 0;
double eliminate_loops_zaxis_cutoff = 2e-2;

/*
 * configuration of vortex injection
 */
//injecting loops at upper z-plane
int loop_injection = 0; //bool, inject or not
double loop_injection_frequency = 1; //injections per second

//line pair injection
int line_injection = 0; //bool, inject or not
int line_injection_n = 1; //how many pairs to inject
double line_injection_frequency; //injections per second
int line_injection_polarized = 0; //whether to inject the vortices in a polarized way, off by default

//Barnes-Hut
int use_BH = 0; //false by default
double BH_resolution = 1e-3;
